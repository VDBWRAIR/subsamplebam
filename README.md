[![Build Status](https://travis-ci.org/necrolyte2/subsamplebam.svg?branch=master)](https://travis-ci.org/necrolyte2/subsamplebam)
[![Coverage Status](https://coveralls.io/repos/necrolyte2/subsamplebam/badge.png?branch=master)](https://coveralls.io/r/necrolyte2/subsamplebam?branch=master)
[![Docs](https://readthedocs.org/projects/subsamplebam/badge/?version=latest)](http://subsamplebam.readthedocs.org/en/latest/)

# subsamplebam

Hopefully will allow you to randomly subsample a bam to get closer to even depth across the genome

Here are some quick comparisons:

## The original sample looks like this

![original](/images/original.png)

## Using samtools to get a subsample

Here we have a bam file called 8457 which has over 5 million reads in it. I picked 1.02 to try and get about 100k of those reads 

```
$> samtools view -hb -s 1.02 8457.bam > smaller.bam; samtools index smaller.bam
```

![samtools](/images/samtools.png)

## Using subsamplebam.py

```
$> samtools view -H 8457.bam
@HD VN:1.3  SO:coordinate
@SQ SN:Den1/GU131895_1/Cambodia/2009/Den1_1 LN:10474
@RG ID:Roche454 SM:8457 CN:None PL:L454
@RG ID:IonTorrent   SM:8457 CN:None PL:IONTORRENT
@RG ID:MiSeq    SM:8457 CN:None PL:ILLUMINA
@RG ID:Sanger   SM:8457 CN:None PL:CAPILLARY
$> subsamplebam 8457.bam Den1/GU131895_1/Cambodia/2009/Den1_1:1-10747 --subsample 10 | samtools view -hSb - | samtools sort - subsampled; samtools index subsampled.bam;
```

![subsamplebam](/images/subsamplebam.png)

## Understanding the output

subsamplebam outputs to both stderr and stdout.
The output to stderr will be the currect samtools region it is working on as well as sometimes you may see a line about Depth being only a certain amount.

Example output:

```
samtools view smaller.bam Den1/GU131895_1/Cambodia/2009/Den1_1:4600-4600
samtools view smaller.bam Den1/GU131895_1/Cambodia/2009/Den1_1:5744-5744
Depth for Den1/GU131895_1/Cambodia/2009/Den1_1:3693-3693 is only 4
samtools view smaller.bam Den1/GU131895_1/Cambodia/2009/Den1_1:3694-3694
samtools view smaller.bam Den1/GU131895_1/Cambodia/2009/Den1_1:4335-4335
samtools view smaller.bam Den1/GU131895_1/Cambodia/2009/Den1_1:5467-5467
```

The lines that begin with samtools are just showing you what position is being evaluated
The lines that begin with Depth are telling you that at that location, the read depth was only X which was lower than what you specified with ``--subsample``

## Installing samtools

If you don't already have samtools installed and in your PATH you can easily install it into this project as follows

```
tests/install_samtools.sh
export PATH=$PATH:$PWD/bin
```

This will only work for your current terminal session unless you modify your PATH in your .bashrc
