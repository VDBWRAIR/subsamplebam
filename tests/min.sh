#!/bin/bash
# run my minimization script and create the plots
cd ..
dir="foo"
mkdir -p $dir
tmp=${dir}/${2}.min.${6}
echo $tmp
echo `pwd`
ngs="/home/AMED/michael.panciera/projects/ngs_mapper/ngs_mapper"
python subsample_mindepth.py $1 $2 $3 $4 $5 $6 | samtools view -hSb - | samtools sort - $tmp; samtools index $tmp.bam;
python $ngs/bam_to_qualdepth.py $tmp.bam > $tmp.json
python  $ngs/graph_qualdepth.py $tmp.json -o $tmp.png
