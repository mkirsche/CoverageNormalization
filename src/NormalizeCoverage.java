/*
 * Takes a sorted SAM file and downsamples the reads in a way that tries to get close to uniform coverage
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;

public class NormalizeCoverage {
	
	// Maximum length allowed - make sure it's longer than genome length
	static int MAX_LEN = 31000;
	
	// The amount of coverage we want to keep everywhere
	static int COV_THRESHOLD = 50;
	
	// Whether or not to sort by quality values
	static boolean QUAL_SORT = false;
	
	static int RAND_SEED = -1;
	
	// Whether or not to use Rampart-style CSV files in place of SAM files
	static boolean INPUT_CSV = false;
	
	// Input and output filenames
	static String fn = "", ofn = "";
	
	// File containing the cumulative coverage information with one byte per position
	static String coverageFn = "";
	
	// Whether or not to perform logging
	static boolean logStats = true;
	
	/*
	 * Prints out usage instructions
	 */
	static void usage()
	{
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
		System.out.println("  covfile            (String) []      - the file containing coverage from other samples");
		System.out.println("  --qual_sort                         - prioritize reads with higher alignment quality");
		System.out.println("  --input_csv                         - expect the input to be a Rampart-formatted CSV file");
		System.out.println("  --no_logging                        - don't produce logging files");

		System.out.println();
	}
	
	/*
	 * Parses command line arguments
	 */
	static void parseArgs(String[] args)
	{
		for(String str : args)
		{
			String s = str;
			int equalsIdx = s.indexOf('=');
			if(equalsIdx == -1)
			{
				if(s.endsWith("qual_sort"))
				{
					QUAL_SORT = true;
				}
				if(s.endsWith("input_csv"))
				{
					INPUT_CSV = true;
				}
				if(s.endsWith("no_logging"))
				{
					logStats = false;
				}
			}
			else
			{
				String key = s.substring(0, equalsIdx).toLowerCase();
				String val = s.substring(1 + equalsIdx);
				
				if(key.equals("input"))
				{
					fn = val;
				}
				else if(key.equals("output"))
				{
					ofn = val;
				}
				else if(key.equals("covfile"))
				{
					coverageFn = val;
				}
				else if(key.equals("coverage_threshold"))
				{
					COV_THRESHOLD = Integer.parseInt(val);
				}
				else if(key.equals("genome_max_len"))
				{
					MAX_LEN = Integer.parseInt(val);
				}
			}
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		// Print help menu for -h or --help
		if(args.length == 0 || args[0].equalsIgnoreCase("-h") || args[0].equalsIgnoreCase("--help"))
		{
			usage();
			System.exit(0);
		}
		else
		{
			parseArgs(args);
		}
		
		// Check that input file actually exists
		if(!new File(fn).exists())
		{
			System.out.println("\nInput file " + fn + " does not exist!\n");
			System.exit(1);
		}
		
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		
		// Get the reference intervals for all reads
		ArrayList<Read> reads = new ArrayList<Read>();
		int lineIdx = 0;
		while(input.hasNext())
		{
			String line = input.nextLine();
			
			// Ignore SAM header lines
			boolean isReadLine = true;
			if(!INPUT_CSV && line.startsWith("@"))
			{
				isReadLine = false;
			}
			
			// Ignore first line in CSV
			if(INPUT_CSV && lineIdx == 0)
			{
				isReadLine = false;
			}
			
			// Increment line number
			lineIdx++;
			
			if(!isReadLine)
			{
				continue;
			}
			
			// Initialize read fields: start and end on reference, read length, and quality score
			int[] startEnd = new int[2];
			int rl = 0;
			double qual = 0.0;
			
			if(INPUT_CSV)
			{
				// CSV format - parse fields
				String[] tokens = line.split(",");
				
				startEnd = new int[] {Integer.parseInt(tokens[6]), Integer.parseInt(tokens[7])};
				rl = Integer.parseInt(tokens[1]);
				qual = 1.0 * Integer.parseInt(tokens[8]) / rl;
			}
			else
			{
				// SAM format - parse fields
				String[] tokens = line.split("\t");

				startEnd = refInterval(line);
				rl = cigarQueryLength(tokens[5]);
				qual = 1.0 * cigarNumMatches(tokens[5]) / rl;
			}
			
			// Add the read to the list	
			reads.add(new Read(reads.size(), startEnd[0], startEnd[1], rl, qual));
		}
		
		int n = reads.size();
		
		// This will be filled with total coverage of each position
		int[] cov = new int[MAX_LEN];
		
		// Add +1 to represent coverage going up at start and -1 to represent coverage down at end
		for(int i = 0; i<n; i++)
		{
			Read r = reads.get(i);
			cov[r.start]++;
			cov[r.end]--;
		}
		
		// Now each element of cov will be coverage(i) - coverage(i-1),
		// so take the cumulative sum to make it actual coverage
		for(int i = 1; i<cov.length; i++)
		{
			cov[i] += cov[i-1];
		}
		
		if(QUAL_SORT)
		{
			Collections.sort(reads, new Comparator<Read>() {

				@Override
				public int compare(Read a, Read b) {
					// TODO Auto-generated method stub
					return Double.compare(b.qual, a.qual);
				}});
		}
		else
		{
			Collections.shuffle(reads);
		}
		
		// The coverage so far of each position by reads we choose to keep
		int[] readCov = new int[MAX_LEN];
		
		// Initialize reader for coverage file
		RandomAccessFile coverageFileReader = null;
		if(coverageFn != null && coverageFn.length() > 0 && new File(coverageFn).exists())
		{
			if(new File(coverageFn).length() > 0)
			{
				coverageFileReader = new RandomAccessFile(coverageFn, "r");
			}
		}
		
		// If we have a coverage file, add the data there to the read coverage array
		if(coverageFileReader != null)
		{
			byte[] buf = new byte[MAX_LEN];
			coverageFileReader.seek(0);
			coverageFileReader.read(buf);
			for(int i = 0; i<MAX_LEN; i++) readCov[i] += buf[i];
		}
		
		// True for reads we want to keep
		boolean[] used = new boolean[n];
		
		// Go through the reads, and if there's some position covered by it that's
		// below coverage threshold, take the read
		for(Read r : reads)
		{
			// Check the minimum coverage across the interval of the read
			int minCov = readCov[r.start];
			for(int i = r.start; i < r.end; i++)
			{
				minCov = Math.min(readCov[i], minCov);
			}
			
			// If the minimum is low enough, take the read and add its coverage
			if(minCov < COV_THRESHOLD)
			{
				used[r.index] = true;
				for(int i = r.start; i < r.end; i++)
				{
					readCov[i]++;
				}
			}
		}
		
		// Calculate some statistics
		
		// Total quality score of all reads
		double totalQual = 0;
		
		// Total quality score of all used reads
		double usedTotalQual = 0;
		
		// Total number of reads
		int totalCount = n;
		
		// Total number of used reads
		int usedCount = 0;
		
		// Total number of bases across all reads
		int totalBases = 0;
		
		// Total number of bases across all used reads
		int usedBases = 0;
		
		// Minimum coverage in full dataset (ignoring first and last 50 bp)
		int minCov = 987654321;
		
		// Minimum coverage among kept reads (ignoring first and last 50 bp)
		int totalMin = 987654321;
		
		// Compute the statistics outlind above
		for(int i = 0; i<n; i++)
		{
			Read r = reads.get(i);
			totalBases += r.end - r.start;
			totalQual += r.qual;
			if(used[r.index])
			{
				usedTotalQual += r.qual;
				usedCount++;
				usedBases += r.end - r.start;
			}
		}
		
		// Check for min coverage outside of first/last 50 bp and see if it's similar to old minimum
		for(int i = 50; i<MAX_LEN-50; i++)
		{
			if(cov[i] > 0)
			{
				minCov = Math.min(minCov, cov[i]);
				totalMin = Math.min(totalMin, readCov[i]);
			}
			
			if(cov[i] >= COV_THRESHOLD && readCov[i] < COV_THRESHOLD)
			{
				System.out.println("Coverage dropped below threshold at position " + i+"; "
						+ "Old coverage="+cov[i]+", New coverage="+readCov[i]);
			}
		}
		
		// Check for any differences in coverage in the 50 bases on each end, which I'm assuming to have lower coverage
		for(int i = 0; i<50; i++)
		{
			if(cov[i] != readCov[i])
			{
				System.out.println("Uneven coverage near ends at position " + i+"; "
						+ "Old coverage =" + cov[i]+", New coverage="+readCov[i]);
			}
			if(cov[MAX_LEN - i - 1] != readCov[MAX_LEN - i - 1])
			{
				System.out.println("Uneven coverage near ends at position " + (MAX_LEN - i - 1)+"; "
						+ "Old coverage =" + cov[MAX_LEN - i - 1]+", New coverage="+readCov[MAX_LEN - i - 1]);
			}
		}
		
		// Output statistics
		System.out.println("Total read count (unfiltered): " + totalCount);
		System.out.println("Downsampled read count: " + usedCount);
		
		System.out.println("Overall average alignment accuracy: " + String.format("%.6f", 1.0 * totalQual / totalCount));
		System.out.println("Downsampled average alignment accuracy: " + String.format("%.6f", 1.0 * usedTotalQual / usedCount));
		System.out.println("Total bases covered (unfiltered): " + totalBases);
		System.out.println("Downsampled bases covered: " + usedBases);
		
		System.out.println("Old min coverage: " + totalMin);
		System.out.println("Downsampled min coverage: " + minCov);
		
		// Go through reads and make file with filtered reads
		
		input = new Scanner(new FileInputStream(new File(fn)));
		
		// Generate output filename
		if(ofn.length() == 0)
		{
			String suff = INPUT_CSV ? ".csv" : ".sam";
			if(fn.endsWith(suff))
			{
				ofn = fn.substring(0, fn.length() - 4) + ".covfiltered" + suff;
			}
			else
			{
				ofn = fn + ".covfiltered" + suff;
			}
		}
		
		// Write out the reads we want to keep
		PrintWriter out = new PrintWriter(new File(ofn));
		int readIndex = 0;
		lineIdx = 0;
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(!INPUT_CSV && line.startsWith("@"))
			{
				out.println(line);
				lineIdx++;
				continue;
			}
			if(INPUT_CSV && lineIdx == 0)
			{
				out.println(line);
				lineIdx++;
				continue;
			}
			lineIdx++;
			if(used[readIndex])
			{
				out.println(line);
			}
			readIndex++;
		}
		
		// Output coverage at each base if a coverage file was provided
		if(coverageFn != null && coverageFn.length() > 0)
		{
			coverageFileReader = new RandomAccessFile(coverageFn, "rw");
			byte[] toWrite = new byte[MAX_LEN];
			for(int i = 0; i<MAX_LEN; i++)
			{
				toWrite[i] = (byte)Math.min(readCov[i], COV_THRESHOLD);
			}
			coverageFileReader.write(toWrite, 0, MAX_LEN);
			coverageFileReader.close();
		}
		
		if(logStats)
		{
			// Print the old and new coverage of each base
			String coverageReadableFn = "coverage.txt";
			PrintWriter coverageReadableOut = new PrintWriter(new File(coverageReadableFn));
			
			for(int i = 0; i<MAX_LEN; i++)
			{
				if(cov[i] > 0)
				{
					coverageReadableOut.println(cov[i]+" "+readCov[i]);
				}
			}
			coverageReadableOut.close();
			
			// Print out the read length in the whole dataset and in the sample
			String allLengthsFn = "lengths_all.txt";
			String sampleLengthsFn = "lengths_sample.txt";
			PrintWriter allLengthsOut = new PrintWriter(new File(allLengthsFn));
			PrintWriter sampleLengthsOut = new PrintWriter(new File(sampleLengthsFn));
			for(int i = 0; i<n; i++)
			{
				allLengthsOut.println(reads.get(i).readLength);
				if(used[reads.get(i).index])
				{
					sampleLengthsOut.println(reads.get(i).readLength);
				}
			}
			allLengthsOut.close();
			sampleLengthsOut.close();
		}
		
		input.close();
		out.close();
	}
	
	/*
	 * Read represented by line number in the SAM file and start/end positions in the reference
	 */
	static class Read implements Comparable<Read>
	{
		int index, start, end;
		int readLength;
		double qual;
		
		Read(int ii, int ss, int ee, int rl, double qq)
		{
			index = ii;
			start = ss;
			end = ee;
			readLength = rl;
			qual = qq;
		}

		@Override
		public int compareTo(Read o) {
			// TODO Auto-generated method stub
			return start - o.start;
		}
	}
	
	/*
	 * Gets the (start, end) in the reference of a SAM-format line
	 */
	static int[] refInterval(String line)
	{
		String[] tokens = line.split("\t");
		int start = Integer.parseInt(tokens[3]);
		
		String cigar = tokens[5];
		int end = start + cigarRefLength(cigar);
		return new int[] {start, end};
	}
	
	/*
	 * Whether or not a CIGAR character consumes the reference
	 */
	static boolean consumesReference(char c)
	{
		return c == 'M' || c == 'D' || c == 'N' || c == '=' || c == 'X';
	}
	
	static int cigarNumMatches(String s)
	{
		int curLen = 0;
		int res = 0;
		for(int i = 0; i<s.length(); i++)
		{
			char c = s.charAt(i);
			if(c >= '0' && c <= '9')
			{
				curLen = curLen * 10 + (c - '0');
			}
			else
			{
				if(c == 'M' || c == '=')
				{
					res += curLen;
				}
				curLen = 0;
			}
		}
		return res;
	}
	
	/*
	 * Whether or not a CIGAR character consumes the reference
	 */
	static boolean consumesQuery(char c)
	{
		return c == 'M' || c == 'I' || c == 'S' || c == '=' || c == 'X';
	}
	
	static int cigarQueryLength(String s)
	{
		int curLen = 0;
		int res = 0;
		for(int i = 0; i<s.length(); i++)
		{
			char c = s.charAt(i);
			if(c >= '0' && c <= '9')
			{
				curLen = curLen * 10 + (c - '0');
			}
			else
			{
				if(consumesQuery(c))
				{
					res += curLen;
				}
				curLen = 0;
			}
		}
		return res;
	}
	
	/*
	 * Gets the length of reference spanned by a cigar string
	 */
	static int cigarRefLength(String s)
	{
		int curLen = 0;
		int res = 0;
		for(int i = 0; i<s.length(); i++)
		{
			char c = s.charAt(i);
			if(c >= '0' && c <= '9')
			{
				curLen = curLen * 10 + (c - '0');
			}
			else
			{
				if(consumesReference(c))
				{
					res += curLen;
				}
				curLen = 0;
			}
		}
		return res;
	}
}
