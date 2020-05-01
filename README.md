# CoverageNormalization

## Description

This script takes a list of alignments in SAM format and removes reads from regions with higher coverage to obtain a smaller sample with more uniform coverage.

It uses a threshold-based algorithm where it orders the reads either randomly (default) or in order of alignment accuracy, and incrementally keeps reads which have at least one base that has not yet reached a user-defined coverage threshold.

This gives the algorithm a few properties:
* For any position in the genome with total coverage less than the threshold, all reads spanning it will be kept
* For any position in the genome with total coverage greater than or equal to the threshold, its coverage among the kept reads will not drop below the threshold.

Slides giving more details about the method as well as preliminary results can be found here: [CoverageNormalization.pdf](https://github.com/mkirsche/CoverageNormalization/blob/master/slides/CoverageNormalization.pdf)

## Compilation

``javac src/*.java``

## Running

```
Usage: java -cp src NormalizeCoverage [args]
  Example: java -cp src NormalizeCoverage input=jhu004.sam

Required args:
  input (String) - a SAM file with the alignments of all of the reads
  
Optional args:
  coverage_threshold (int)    [50]    - the coverage to require at each base (if original coverage is high enough)
  genome_max_len     (int)    [31000] - an upper bound on the genome length
  output             (String) []      - the file to write downsampled reads to
  covfile            (String) []      - the file containing coverage from other samples
  --qual_sort                         - prioritize reads with higher alignment quality
  --input_csv                         - expect the input to be a Rampart-formatted CSV file
  --no_logging                        - don't produce logging files
  --even_strand                       - tries to get even coverage between the strands when possible
  ```
  
## Other Scripts

### Post-processing and plotting
* *plot_coverage.py* - Generate bar plots and histograms of the coverage of each position before and after normalization
* *plotReadLengths.py* - Plots histograms of read lengths before and after normalization
* *PlotStrandBias.py* - Plots histograms of strand bias on full dataset and two samples

### Pipelines/Utilities
* *run_medaka_sample.sh* - Runs medaka for consensus and variant calling on a particular BAM file
* *src/ParseResults.java* - Parses the output of a number of medaka runs and counts the frequency of each variant
* *testdownsampling.sh* - Produces 10 samples with a given threshold and runs medaka on each, counting the frequency of each variant
* *variant_heatmap.py* - Plots heatmaps of the number of trials where different SNPs were called with medaka
* *run_full_test.sh* - Tests 10 samples for each of a number of coverage thresholds and generates heatmaps for all variant calls and for homozygous variant calls
* *strandbiastest.sh* - Runs normalization with and without strand adjustment and tests difference in resulting strand bias
* *src/Vcf2Fasta.java* - Takes a reference and a VCF, and produces a consensus genome sequence with all of the given variants integrated

