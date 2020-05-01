javac src/*.java
java -cp src NormalizeCoverage input=jhu004.sam --qual_sort --even_strand > newstats.txt
mv strandbiassample.txt strandbiassample_new.txt
java -cp src NormalizeCoverage input=jhu004.sam --qual_sort > oldstats.txt
mv strandbiassample.txt strandbiassample_old.txt
python PlotStrandBias.py strandbiasfull.txt strandbiassample_old.txt strandbiassample_new.txt sb.png
