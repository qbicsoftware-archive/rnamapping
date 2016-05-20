import os
import sys
import subprocess
import json
from os.path import join as pjoin
from os.path import exists as pexists
import glob
import numpy as np
import pandas as pd

configfile: "config.json"
workdir: config["var"]
SNAKEDIR = config['src']

try:
    VERSION = subprocess.check_output(
        ['git', 'describe', '--tags', '--always', '--dirty'],
        cwd=SNAKEDIR
    ).decode().strip()
except subprocess.CalledProcessError:
    VERSION = 'unknown'

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
ETC = config['etc']

def data(path):
    return os.path.join(DATA, path)

def ref(path):
    return os.path.join(REF, path)

def log(path):
    return os.path.join(LOGS, path)

def result(path):
    return os.path.join(RESULT, path)

def etc(path):
    return os.path.join(ETC, path)


try:
    with open(etc("params.json")) as f:
        parameters = json.load(f)
except OSError as e:
    print("Could not read parameter file: " + str(e), file=sys.stderr)
    sys.exit(1)
except ValueError as e:
    print("Invalid parameter file: " + str(e), file=sys.stderr)
    sys.exit(1)

default_params = {
    "stranded": 'no',
    "overlap_mode": 'union',
    "gff_attribute": 'gene_id',
    "feature_type": 'exon',
}
default_params.update(parameters)
parameters = default_params

for key in ['gtf', 'stranded', 'overlap_mode', 'indexed_genome',
            'gff_attribute', 'feature_type']:
    if key not in parameters:
        print("Missing parameter %s in etc/params.json" % key, file=sys.stderr)
        exit(1)

parameters['indexed_genome'] = ref(parameters['indexed_genome'])
parameters['gtf'] = ref(parameters['gtf'])

indexed_genome = parameters["indexed_genome"]
if not os.path.exists(indexed_genome + '.fa'):
    raise ValueError("Could not find indexed genome file %s" % indexed_genome)


INPUT_FILES = []
for name in os.listdir(DATA):
    if name.lower().endswith('.sha256sum'):
        continue
    if name.lower().endswith('.fastq'):
        if not name.endswith('.fastq'):
            print("Extension fastq is case sensitive.", file=sys.stderr)
            exit(1)
        INPUT_FILES.append(os.path.basename(name)[:-6])
    elif name.lower().endswith('.fastq.gz'):
        if not name.endswith('.fastq.gz'):
            print("Extension fastq is case sensitive.", file=sys.stderr)
            exit(1)
        INPUT_FILES.append(os.path.basename(name)[:-len('.fastq.gz')])
    elif name.lower().endswith('.bam'):
        if not name.endswith('.bam'):
            print("Extension bam is case sensitive.", file=sys.stderr)
            exit(1)
    elif name.endswith('.bam.bai'):
        continue
    else:
        print("Unknown data file: %s" % name, file=sys.stderr)
        exit(1)

if len(set(INPUT_FILES)) != len(INPUT_FILES):
    print("Some input file names are not unique")
    exit(1)


DESIGN = pd.read_csv(etc('GROUPS'), sep='\t')

RUNS = {}
for group, df in DESIGN.groupby('group'):
    if '_R1' in df['file'].iloc[0] or '_R2' in df['file'].iloc[0]:
        first_id = '_R1'
        second_id = '_R2'
    elif '_1' in df['file'].iloc[0] or '_2' in df['file'].iloc[0]:
        first_id = '_1'
        second_id = '_2'
    else:
        raise ValueError("Unknown paired end sequencing identifiers")

    first = df[df['file'].str.contains(first_id)]
    first = first.sort_values(by='file').reset_index(drop=True)
    second = df[df['file'].str.contains(second_id)]
    second = second.sort_values(by='file').reset_index(drop=True)
    assert(len(df) == len(first) + len(second))
    first = list(sorted(first['file'].values))
    second = list(sorted(second['file'].values))
    prefix = [name[:name.find(first_id)] for name in first]
    for prefix, file1, file2 in zip(prefix, first, second):
        assert (group, prefix) not in RUNS
        assert prefix and group
        assert '___' not in group and '___' not in prefix
        RUNS[(group, prefix)] = (file1, file2)

OUTPUT_FILES = []
OUTPUT_FILES.extend(expand(result("fastqc/{name}"), name=INPUT_FILES))
#OUTPUT_FILES.extend(expand(result("FastQCcut_{name}.zip"), name=INPUT_FILES))
#OUTPUT_FILES.extend(expand(result("/HTSeqCounts_{name}.txt"), name=INPUT_FILES))


#OUTPUT_FILES.extend(expand("Summary/NumReads/Original/{name}.txt", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.extend(expand("Summary/NumReads/PreFilter/{name}.txt", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.extend(expand("{result}/FastQC_{name}.zip", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.extend(expand("{result}/FastQCcut_{name}.zip", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.extend(expand("Summary/NumReads/CutAdaptMerge/{name}.txt", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.extend(expand("{result}/HTSeqCounts_{name}.txt", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.extend(expand("TopHat2/{name}/accepted_hits.bai", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.extend(expand("Summary/MappingStats/{name}.txt", name=INPUT_FILES, result=RESULT))
#OUTPUT_FILES.append("checksums.ok")
#OUTPUT_FILES.append(result('all_counts.csv'))


rule all:
    input: result("all_counts.csv"), OUTPUT_FILES, "checksums.ok"

rule checksums:
    output: "checksums.ok"
    threads: 1
    run:
        out = os.path.abspath(str(output))
        cksums = glob.glob(data("*.sha256sum"))
        if cksums:
            shell("cd %s; "
                  "sha256sum -c *.sha256sum && "
                  "touch %s" % (data('.'), out))
        else:
            shell("touch %s" % out)

rule LinkUncompressed:
    input: data("{name}.fastq")
    output: "fastq/{name}.fastq"
    shell: "ln -s {input} {output}"

rule Uncompress:
    input: data("{name}.fastq.gz")
    output: "fastq/{name}.fastq"
    shell: "zcat {input} > {output}"

rule PreFilterReads:
    input: "fastq/{name}.fastq"
    output: "PreFilterReads/{name}.fastq"
    run:
        with open(str(input)) as infile, open(str(output), 'w') as outfile:
            num_ok = 0
            num_filtered = 0
            line = infile.readline()
            while line:
                assert line.startswith('@')
                body = ''.join(infile.readline() for _ in range(3))
                if all(':Y:' not in part for part in line.split(' ')[1:]):
                    num_ok += 1
                    outfile.write(line)
                    outfile.write(body)
                else:
                    num_filtered += 1
                line = infile.readline()
        num_all = num_ok + num_filtered
        if num_filtered / num_all > .1:
            raise ValueError("More than 10% of reads were filtered in "
                "PreFilterReads. This probably indicates a bug.")

rule fastqc:
    input: "PreFilterReads/{name}.fastq"
    output: result("fastqc/{name}")
    threads: 1
    run:
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        shell("fastqc {input} -o {output}")

rule FastQCcut:
    input: "CutAdaptMerge/{name}.fastq"
    output: result("fastqccut/{name}")
    threads: 1
    run:
        try:
            os.mkdir(str(output))
        except Exception:
            pass
        shell("fastqc {input} -o {output}")

rule Overrepresented:
    input: result("fastqc/{name}")
    output: "Overrepresented/{name}.txt"
    run:
        f = open(str(input) + "/" + wildcards['name'] + "_fastqc/fastqc_data.txt")
        out = open(str(output), "w")
        sw = False
        for line in f:
            if line.startswith(">>Overrepresented"):
                sw = True
            if not line.startswith(">>END_MODULE") and sw:
                out.write(line)
            else:
                sw = False
        f.close()
        out.close()

rule OverrepTxtFasta:
    input: "Overrepresented/{name}.txt"
    output: "Overrepresented/{name}.fasta"
    run:
        f = open(str(input))
        out = open(str(output), "w")
        for line in f:
            if not (line.startswith("#") or line.startswith(">>")):
                 tokens1 = line.split("\t")
                 print(tokens1)
                 fseq = tokens1[0]
                 fheader = '>' + tokens1[3]
                 out.write(fheader)#+'\n')#it just happens that is the last column
                 out.write(fseq+'\n')
        f.close()
        out.close()

rule MergeAdapters:
    input: expand("Overrepresented/{name}.fasta", name=INPUT_FILES)
    output: "MergeAdapters/merged.fasta"
    shell: "cat {input} > {output}"

rule CutAdapt:
    input: 
        adapter="MergeAdapters/merged.fasta",
        fastq=lambda w: ['PreFilterReads/' + p for p in RUNS[(w['group'], w['prefix'])]]
        #fastq=lambda w: [data(p) for p in RUNS[(w['group'], w['prefix'])]]
    output: 
        L="CutAdaptMerge/{group}___{prefix}_R1.fastq",
        R="CutAdaptMerge/{group}___{prefix}_R2.fastq"
    run:
        with open(input['adapter']) as f:
            skip = not bool(f.read(1))
        if skip:
            os.symlink(os.path.abspath(str(input['adapter'])), str(output))
        else:
            shell('cutadapt --discard-trimmed -a file:' + input['adapter'] + ' -o ' + output['L'] + ' -p ' + output['R'] + ' ' + ' '.join(input['fastq']))

rule TopHat2:
    #input: lambda w: ['CutAdaptMerge/' + p for p in RUNS[(w['group'], w['prefix'])]]
    input: 
        L="CutAdaptMerge/{group}___{prefix}_R1.fastq",
        R="CutAdaptMerge/{group}___{prefix}_R2.fastq"
    output: "TopHat2/{group}___{prefix}"
    run:
        gtf = parameters['gtf']
        genome = parameters['indexed_genome']
        shell('tophat --no-coverage-search -o {output}'+ ' -p 2 -G %s %s %s %s'
              % (gtf, genome, input['L'], input['R']))

rule HTSeqCounts:
    input: "TopHat2/{group}___{prefix}"
    output: result("HTSeqCounts_{group}___{prefix}.txt")
    run:
        sam_command = "samtools view {input}/accepted_hits.bam"
        htseq = ("htseq-count -i {gff_attribute} -t {feature_type} "
                 "-m {overlap_mode} -s {stranded} - {gtf}").format(**parameters)
        shell("%s | %s > {output}" % (sam_command, htseq))

rule CombineCounts:
    input: expand(result("HTSeqCounts_{group}___{prefix}.txt").format(group=key[0], prefix=key[1]) for key in RUNS)
    output: result("all_counts.csv")
    run:
        import re

        pattern = "HTSeqCounts_([0-9a-zA-Z_\- ]*).txt"
        names = [re.search(pattern, str(name)).groups()[0] for name in input]
        data = {}
        for name, file in zip(names, input):
            file = str(file)
            data[name] = pd.Series.from_csv(file, sep='\t')
        df = pd.DataFrame(data)
        df.index.name = parameters['gff_attribute']
        df.to_csv(str(output))

rule IndexBAM:
    input: "TopHat2/{name}"
    output: "TopHat2/{name}/accepted_hits.bai"
    shell: "samtools index {input}/accepted_hits.bam && mv -f {input}/accepted_hits.bam.bai  {output}"

rule CpAlignSummary:
    input: "TopHat2/{name}"
    output: "Summary/MappingStats/{name}.txt"
    shell: "cp {input}/align_summary.txt {output}"

rule PerBaseCoverage:
    input: "TopHat2/{name}"
    output: "Statistics/{name}_PerBaseCoverage.txt"
    shell: "samtools depth {input}/accepted_hits.bam > {output}"

rule Numreads:
    input: "PreFilterReads/{name}.fastq"
    output: "Summary/NumReads/PreFilter/{name}.txt"
    shell: '''dc -e "$(wc -l {input} | cut -f1 -d' ') 4 / p" > {output}'''

rule NumreadsCut:
    input: "CutAdaptMerge/{name}.fastq"
    output: "Summary/NumReads/CutAdaptMerge/{name}.txt"
    shell: '''dc -e "$(wc -l {input} | cut -f1 -d' ') 4 / p" > {output}'''

rule NumreadsOrig:
    input: "fastq/{name}.fastq"
    output: "Summary/NumReads/Original/{name}.txt"
    shell: '''dc -e "$(wc -l {input} | cut -f1 -d' ') 4 / p" > {output}'''
