RNA-Seq workflow for paired-end data. The workflow can be downloaded and run on a cluster environment using the tool qproject provided on github: https://github.com/qbicsoftware/qproject

The workflow uses a module system to load the required software. Be sure to check the `jobscript.sh` file to see which software modules are required. The modules are loaded automatically when using `qproject run` to start the workflow, otherwise they have to be loaded manually.

One should use qproject to download the files. This also creates all folders necessary for the workflow.

```
qproject create -t . -w github:qbicsoftware/rnaseq
```

Be sure to add a config file "params.json" in `etc` which should look like this:

```json
{
    "gtf": "path/to/gtf",
    "indexed_genome": "path/to/genome/basename",
    "stranded": "yes",
    "overlap_mode": "union",
    "feature_type": "exon",
    "gff_attribute": "gene_id",
    "normalize_counts": "deseq2"
}
```

where `indexed_genome` and `gtf` are paths relative to `ref`.

`indexed_genome` is the basename of a bowtie2 index.
`gtf` is the .gtf file of the reference genome

The parameters `stranded`, `overlap_mode`, `feature_type` and `gff_attribute`
are explained in the htseq documentation.

Members of QBiC can download the data for analysis using `qpostman`: https://github.com/qbicsoftware/postman-cli

It should be installed on the computing stations and can be loaded with:

```
module load qbic/qpostman/0.1.2.3
```

To download the data navigate to the data folder and either provide a QBiC ID
```
java -jar qpostman.jar -i <QBiC ID> -u <your_qbic_username>
```

or a file containing the project QBiC IDs:

```
postman-0.1.2.3 -f sample_IDs.csv -u user-name
```

If you're not using `qpostman` just put the relevant files in the data folder (formats supported: `.fastq`, `.fastq.gz`).

Note that the paired end data should have a naming similar to C13R014c02_05_S6_L008_001_R1.fastq, C13R014c02_05_S6_L008_001_R2.fastq containing the identifier R1 and R2 at the end of the file name.

Next one needs to create the `GROUPS` file, which is a file with each group (e.g. C13R014c02 for the previous example) and the corresponding file names (C13R014c02_05_S6_L008_001_R1.fastq and C13R014c02_05_S6_L008_001_R2.fastq).
This file can be created automatically using the shell script `add_GROUPS_file.sh`. For this navigate to the src folder and use:

```
bash add_GROUPS_file.sh
```

After running the shell script the `GROUPS` file should be present in the `etc` folder.


To run the workflow navigate to the `src` folder.
Using `snakemake -n` one can display the operations the workflow will perform.
Using the `--dag` parameter and piping it to `dot` one can create a .pdf version of the directed acyclic graph used by snakemake to inspect the behavious of the workflow on a local machine.

```
cd src/
snakemake -n
snakemake --dag | dot -Tpdf > dag.pdf
```

For running:

```
qproject run -t ..
```

While running one can inspect the log files (e.g. in a different screen session) for the progress and errors generated by the workflow:

```
cd logs/
tail snake.err -f
```

And to check the jobs on the computing cluster one can use `qstat`.

Alternatively to using `qproject run` one could use `snakemake -j` to run the workflow, but then be sure to check the `jobscript.sh` to load the required modules manually and also note that this would also not use `qsub` to submit the jobs.
