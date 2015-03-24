#!/bin/bash
# run my minimization script and create the plots
cd ..
mkdir minz
tmp="minz/recent.947" 
ngs="/home/AMED/michael.panciera/projects/ngs_mapper/ngs_mapper"
python subsample_mindepth.py $1 $2 $3 | samtools view -hSb - | samtools sort - $tmp; samtools index $tmp.bam;
python $ngs/bam_to_qualdepth.py $tmp.bam > $tmp.json
python  $ngs/graph_qualdepth.py $tmp.json -o $tmp.png
