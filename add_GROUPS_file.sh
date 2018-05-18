#go to ws and there to data/ folder

cd ../data

for i in *
do 
    echo $i | awk -F'[_]' '{print $1}' | tr '\n' '\t' >> GROUPS
    echo $i >> GROUPS
done

sed -i '1s/^/group\tfile\n/' GROUPS
mv GROUPS ../etc

