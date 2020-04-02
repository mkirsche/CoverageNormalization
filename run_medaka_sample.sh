#!/bin/bash

if [ "$#" -ne 3 ]
then
  echo "run_medaka_sample.sh bamfile consensus vcfout"
  exit
fi

BAMFILE=$1
CONSENSUS=$2
REF=/home/mkirsche/bin/jhu_genomics_workshop/src/rampart_jhuapl/primer_schemes/nCoV-2019/V1/nCoV-2019.reference.fasta
VCFOUT=$3

if [ ! -r $CONSENSUS ]
then
  medaka consensus $BAMFILE $CONSENSUS
fi
medaka snp $REF $CONSENSUS $VCFOUT
