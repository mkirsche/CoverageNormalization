#!/bin/bash

javac src/*.java
OUTDIR=output_sorted
#./testdownsampling.sh 30 $OUTDIR
#./testdownsampling.sh 50 $OUTDIR
#./testdownsampling.sh 100 $OUTDIR
#./testdownsampling.sh 200 $OUTDIR
#./testdownsampling.sh 500 $OUTDIR
./testdownsampling.sh 1000 $OUTDIR
./testdownsampling.sh 10000 $OUTDIR
python variant_heatmap.py snps_sorted.png
python variant_heatmap.py snps_homo_sorted.png --homozygous
