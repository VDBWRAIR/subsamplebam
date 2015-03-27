#!/bin/bash
# run the subsamplebam solution
cd ..
mkdir random
ngs="/home/AMED/michael.panciera/projects/ngs_mapper/ngs_mapper" 
bamfile="../780/780.bam"
refid=`samtools view -H ../780/780.bam | grep @SQ | head -1 | grep -oP "SN:\K([^\t]+)"`
seqlength=`samtools view -H ../780/780.bam | grep @SQ | head -1 | grep -oP "LN:\K([0-9]+)"`
tmp="random/780.subsampled"

echo ${refid}:1-${seqlength} Trying

subsamplebam $bamfile ${refid}:1-${seqlength} --subsample 10 | samtools view -hSb - | samtools sort - $tmp; samtools index ${tmp}.bam;
python $ngs/bam_to_qualdepth.py $tmp.bam > $tmp.json
python  $ngs/graph_qualdepth.py $tmp.json -o $tmp.png



#samtools view -hSb - | samtools sort - $tmp; samtools index ${tmp}.bam;
