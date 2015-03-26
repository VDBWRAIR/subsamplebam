bamfile=$1
depth=$2
orphans=$3
refs=`samtools view -H $bamfile | grep @SQ |  grep -oP "SN:\K([^\t]+)"`
outdir=outputs
tmpdir=/tmp/subsampletmp
mkdir -p $tmpdir
mkdir -p $outdir
outfile=${bamfile//\//_}.minimized.$depth
ngs="/home/AMED/michael.panciera/projects/ngs_mapper/ngs_mapper" 
compiled=$outdir/compiled.${outfile}.bam;  

let i=0
for ref  in $refs; 
do 
    let i++;
    echo $1;
    out=${tmpdir}/${i}; #${outifile}
    python subsample_mindepth.py $bamfile $ref --subsample $depth $orphans > $out.sam;
    samtools view -hSb $out.sam  > $out.bam; #| samtools sort - $out; samtools index $out.bam;
    ref=${ref//[\/\|]/_}; 
    echo "$ref saved to $out.bam";
done; 

if [ $i -ge 2 ];
then  
    echo "merging files."
    samtools merge -f  $compiled  $tmpdir/[0-9]*.bam ; 
else
    mv $tmpdir/$i.bam $compiled
fi


echo "bam files compiled in $outdir/compiled.${outfile}.bam.  Now compiling and plotting";
samtools sort $compiled $compiled;
samtools index  $compiled; 
python $ngs/bam_to_qualdepth.py $compiled > $compiled.json
python  $ngs/graph_qualdepth.py $compiled.json -o $compiled.png 
rm $tmpdir/[0-9]*.[sb]am;
