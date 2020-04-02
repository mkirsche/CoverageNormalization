#!/bin/bash

if [ "$#" -ne 2 ]
then
  echo "testdownsampling.sh coverage outdir"
  exit
fi

COVERAGE=$1
OUTDIR=$2
mkdir $OUTDIR
echo 'Testing downsampling with coverage: '$COVERAGE
echo 'Output directory: '$OUTDIR
for i in {1..10}
do
  echo $i
  java -cp src NormalizeCoverage coverage_threshold=$COVERAGE input=jhu004.sam --qual_sort
  mv jhu004.covfiltered.sam $OUTDIR/jhu004.covfiltered_$i"_"$COVERAGE".sam"
  samtools view -hb $OUTDIR/jhu004.covfiltered_$i"_"$COVERAGE".sam" > $OUTDIR/jhu004.covfiltered_$i.bam
  samtools index $OUTDIR/jhu004.covfiltered_$i.bam
  if [ -r $OUTDIR/consensus_$i.hdf ]
  then
    rm $OUTDIR/consensus_$i.hdf
  fi
  ./run_medaka_sample.sh /home/mkirsche/eclipse-workspace/Covid/$OUTDIR/jhu004.covfiltered_$i.bam /home/mkirsche/eclipse-workspace/Covid/$OUTDIR/consensus_$i.hdf /home/mkirsche/eclipse-workspace/Covid/$OUTDIR/medaka_sample_$i"_"$COVERAGE.vcf
done

for i in `ls /home/mkirsche/eclipse-workspace/Covid/$OUTDIR/medaka_sample_*$COVERAGE.vcf`; do echo $i; cat $i | grep -v '#'; done > /home/mkirsche/eclipse-workspace/Covid/$OUTDIR/vcfs_$COVERAGE.txt

java -cp src ParseResults $COVERAGE $OUTDIR
java -cp src ParseResults $COVERAGE $OUTDIR --homozygous
