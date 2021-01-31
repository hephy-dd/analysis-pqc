
cd VPX33234

for dir in `ls -d HPK_*`; do
    for file in `ls -d $dir/*.txt`; do
        echo $file
        python ../../analysis-pqc/scripts/txt2json.py $file > $file.json
    done
done
