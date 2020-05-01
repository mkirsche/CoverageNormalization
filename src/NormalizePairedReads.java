/*
 * Normalizes read pairs based on their alignments, always keeping both or neither for aligned pairs
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Scanner;

public class NormalizePairedReads {

	// Maximum length allowed - make sure it's longer than genome length
	static int MAX_LEN = 31000;

	// The amount of coverage we want to keep everywhere
	static int COV_THRESHOLD = 50;

	// Whether or not to sort by quality values
	static boolean QUAL_SORT = false;

	static int RAND_SEED = -1;

	// Input and output filenames
	static String fn = "", ofn = "";

	/*
	 * Prints out usage instructions
	 */
	static void usage() {
		System.out.println("Usage: java -cp src NormalizeCoverage [args]");
		System.out.println("  Example: java -cp src NormalizeCoverage input=jhu004.sam");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  input (String) - a SAM file with the alignments of all of the reads");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  coverage_threshold (int)    [50]    - the coverage to require at each base (if original coverage is high enough)");
		System.out.println("  genome_max_len     (int)    [31000] - an upper bound on the genome length");
		System.out.println("  output             (String) []      - the file to write downsampled reads to");
		System.out.println("  --qual_sort                         - prioritize reads with higher alignment quality");
		System.out.println();
	}

	/*
	 * Parses command line arguments
	 */
	static void parseArgs(String[] args) {
		for (String str : args) {
			String s = str;
			int equalsIdx = s.indexOf('=');
			if (equalsIdx == -1) {
				if (s.endsWith("qual_sort")) {
					QUAL_SORT = true;
				}
			} else {
				String key = s.substring(0, equalsIdx).toLowerCase();
				String val = s.substring(1 + equalsIdx);

				if (key.equals("input")) {
					fn = val;
				} else if (key.equals("output")) {
					ofn = val;
				} else if (key.equals("coverage_threshold")) {
					COV_THRESHOLD = Integer.parseInt(val);
				} else if (key.equals("genome_max_len")) {
					MAX_LEN = Integer.parseInt(val);
				}
			}
		}
	}

	public static void main(String[] args) throws Exception {
		// Print help menu for -h or --help
		if (args.length == 0 || args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("--help")) {
			usage();
			System.exit(0);
		} else {
			parseArgs(args);
		}

		// Check that input file actually exists
		if (!new File(fn).exists()) {
			System.out.println("\nInput file " + fn + " does not exist!\n");
			System.exit(1);
		}

		Scanner input = new Scanner(new FileInputStream(new File(fn)));

		// Get the reference intervals for all reads
		HashMap<String, ReadPair> readMap = new HashMap<String, ReadPair>();
		
		int lineCount = 0;
		while (input.hasNext()) {
			String line = input.nextLine();

			// Ignore SAM header lines
			boolean isReadLine = true;
			if (line.startsWith("@")) {
				isReadLine = false;
			}

			if (!isReadLine) {
				continue;
			}
			
			lineCount++;
			
			// SAM format - parse fields
			String[] tokens = line.split("\t");
			
			int flag = Integer.parseInt(tokens[1]);
			if((flag & 4) > 0 || (flag & 256) > 0 || (flag & 2048) > 0)
			{
				// Unmapped or secondary or supplementary alignment - ignore
				continue;
			}

			// Initialize read fields: start and end on reference, read length, and quality
			// score
			int[] startEnd = new int[2];
			int rl = 0;
			double qual = 0.0;
			
			startEnd = NormalizeCoverage.refInterval(line);
			rl = NormalizeCoverage.cigarQueryLength(tokens[5]);
			qual = 1.0 * NormalizeCoverage.cigarNumMatches(tokens[5]) / rl;

			// Add the read to the list
			String readName = tokens[0];
			if(!readMap.containsKey(readName))
			{
				readMap.put(readName, new ReadPair());
			}
			readMap.get(readName).addRead(lineCount-1, startEnd[0], startEnd[1], rl, qual);
		}
		
		ArrayList<ReadPair> reads = new ArrayList<ReadPair>();
		for(String s : readMap.keySet())
		{
			reads.add(readMap.get(s));
		}

		int n = reads.size();

		// This will be filled with total coverage of each position
		int[] cov = new int[MAX_LEN];

		// Add +1 to represent coverage going up at start and -1 to represent coverage
		// down at end
		for (ReadPair r : reads)
		{
			for(int i = 0; i<r.starts.size(); i++)
			{
				cov[r.starts.get(i)]++;
				cov[r.ends.get(i)]--;
			}
		}

		// Now each element of cov will be coverage(i) - coverage(i-1),
		// so take the cumulative sum to make it actual coverage
		for (int i = 1; i < cov.length; i++) {
			cov[i] += cov[i - 1];
		}

		if (QUAL_SORT) {
			Collections.sort(reads, new Comparator<ReadPair>() {

				@Override
				public int compare(ReadPair a, ReadPair b) {
					double aQual = 0.0;
					for(double q : a.quals) aQual += q;
					double bQual = 0.0;
					for(double q : b.quals) bQual += q;
					return Double.compare(bQual, aQual);
				}
			});
		} else {
			Collections.shuffle(reads);
		}

		// The coverage so far of each position by reads we choose to keep
		int[] readCov = new int[MAX_LEN];

		// True for reads we want to keep
		boolean[] used = new boolean[lineCount];

		// Go through the reads, and if there's some position covered by it that's
		// below coverage threshold, take the read
		for (ReadPair r : reads) {
			// Whether or not we want this read
			boolean wantRead = false;
			for(int i = 0; i<r.starts.size(); i++)
			{
				for(int j = r.starts.get(i); j < r.ends.get(i); j++)
				{
					wantRead |= readCov[j] < COV_THRESHOLD;
				}
			}

			// If the minimum is low enough, take the read and add its coverage
			if (wantRead)
			{
				for(int i = 0; i<r.starts.size(); i++)
				{
					used[r.indices.get(i)] = true;
					for (int j = r.starts.get(i); j < r.ends.get(i); j++)
					{
						readCov[j]++;
					}
				}
			}
		}

		// Calculate some statistics

		// Total quality score of all reads
		double totalQual = 0;

		// Total quality score of all used reads
		double usedTotalQual = 0;

		// Total number of reads
		int totalCount = 0;

		// Total number of used reads
		int usedCount = 0;

		// Total number of bases across all reads
		int totalBases = 0;

		// Total number of bases across all used reads
		int usedBases = 0;

		// Minimum coverage in full dataset (ignoring first and last 50 bp)
		int fullMinCov = 987654321;

		// Minimum coverage among kept reads (ignoring first and last 50 bp)
		int sampleMinCov = 987654321;

		// Compute the statistics outlined above
		for (int i = 0; i < n; i++) {
			ReadPair r = reads.get(i);
			for(int j = 0; j<r.starts.size(); j++)
			{
				totalBases += r.ends.get(j) - r.starts.get(j);
				totalQual += r.quals.get(j);
				totalCount++;
				if (used[r.indices.get(j)])
				{
					usedTotalQual += r.quals.get(j);
					usedCount++;
					usedBases += r.ends.get(j) - r.starts.get(j);
				}
			}
		}

		// Check for min coverage outside of first/last 50 bp and see if it's similar to
		// old minimum
		for (int i = 50; i < MAX_LEN - 50; i++) {
			if (cov[i] > 0) {
				fullMinCov = Math.min(fullMinCov, cov[i]);
				sampleMinCov = Math.min(sampleMinCov, readCov[i]);
			}

			if (cov[i] >= COV_THRESHOLD && readCov[i] < COV_THRESHOLD) {
				System.out.println("Coverage dropped below threshold at position " + i + "; " + "Old coverage=" + cov[i]
						+ ", New coverage=" + readCov[i]);
			}

		}

		// Output statistics
		System.out.println("Total read count (unfiltered): " + totalCount);
		System.out.println("Downsampled read count: " + usedCount);

		System.out
				.println("Overall average alignment accuracy: " + String.format("%.6f", 1.0 * totalQual / totalCount));
		System.out.println(
				"Downsampled average alignment accuracy: " + String.format("%.6f", 1.0 * usedTotalQual / usedCount));
		System.out.println("Total bases covered (unfiltered): " + totalBases);
		System.out.println("Downsampled bases covered: " + usedBases);

		System.out.println("Old min coverage: " + fullMinCov);
		System.out.println("Downsampled min coverage: " + sampleMinCov);

		// Go through reads and make file with filtered reads

		input = new Scanner(new FileInputStream(new File(fn)));

		// Generate output filename
		if (ofn.length() == 0) {
			String suff = ".sam";
			if (fn.endsWith(suff)) {
				ofn = fn.substring(0, fn.length() - 4) + ".covfiltered" + suff;
			} else {
				ofn = fn + ".covfiltered" + suff;
			}
		}

		// Write out the reads we want to keep
		PrintWriter out = new PrintWriter(new File(ofn));
		int readIndex = 0;
		while (input.hasNext()) {
			String line = input.nextLine();
			if (line.startsWith("@")) {
				out.println(line);
				continue;
			}
			if (used[readIndex]) {
				out.println(line);
			}
			readIndex++;
		}

		input.close();
		out.close();
	}

	/*
	 * Read represented by line number in the SAM file and start/end positions in
	 * the reference
	 */
	static class ReadPair implements Comparable<ReadPair> {
		ArrayList<Integer> indices;
		ArrayList<Integer> starts;
		ArrayList<Integer> ends;
		ArrayList<Double> quals;
		ArrayList<Integer> readLengths;
		
		ReadPair()
		{
			indices = new ArrayList<Integer>();
			starts = new ArrayList<Integer>();
			ends = new ArrayList<Integer>();
			quals = new ArrayList<Double>();
			readLengths = new ArrayList<Integer>();
		}

		void addRead(int ii, int ss, int ee, int rl, double qq) {
			indices.add(ii);
			starts.add(ss);
			ends.add(ee);
			readLengths.add(rl);
			quals.add(qq);
		}

		@Override
		public int compareTo(ReadPair o) {
			return starts.get(0) - o.starts.get(0);
		}
	}
}
