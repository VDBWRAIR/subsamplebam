#!/bin/bash
set -ev

version=1.1
git clone https://github.com/samtools/htslib.git
git clone https://github.com/samtools/samtools.git
cd samtools
git checkout $version
sed -i "s|prefix\s\+=\s\+/usr/local|prefix = $PWD|" Makefile
make
make install
