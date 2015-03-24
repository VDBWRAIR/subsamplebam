#!/bin/bash

tmp="tmp/minimized.780.depth200"
ngs="/home/AMED/michael.panciera/projects/ngs_mapper/ngs_mapper"
python subsample_mindepth.py $1 $2 | samtools view -hSb - | samtools sort - $tmp; samtools index $tmp.bam;
python $ngs/bam_to_qualdepth.py $tmp.bam > $tmp.json
python  $ngs/graph_qualdepth.py $tmp.json -o $tmp.png
