/*
 * A utility for splicing small variants into a reference genome
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeSet;

public class Vcf2Fasta 
{
	
	static String refFn = "", vcfFn = "", outFn = "";
	
	static int charsPerLine = 60;
	
	static boolean ignoreHeterozygous = true;
	
	/*
	 * Prints out usage instructions
	 */
	static void usage()
	{
		System.out.println("Usage: java -cp src Vcf2Fasta [args]");
		System.out.println("  Example: java -cp src Vcf2Fasta ref=genome.fa vcf=vars.vcf output=genome_new.fa");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  ref    (String) - a FASTA file containing the reference genome");
		System.out.println("  vcf    (String) - a VCF file containing a list of variants (heterozygous variants will be ignored)");
		System.out.println("  output (String) - the FASTA file that will be output with the variants integrated");
		System.out.println("Optional args:");
		System.out.println("  --keep_heterozygous - also splice in variants which are heterozygous (as if they were homozygous)");
		
		System.out.println();
	}
	
	/*
	 * Parses command line arguments
	 */
	static void parseArgs(String[] args) throws Exception
	{
		for(String str : args)
		{
			String s = str;
			int equalsIdx = s.indexOf('=');
			if(equalsIdx == -1)
			{
				if(s.toLowerCase().endsWith("keep_heterozygous"))
				{
					ignoreHeterozygous = false;
				}
			}
			else
			{
				String key = s.substring(0, equalsIdx).toLowerCase();
				String val = s.substring(1 + equalsIdx);
				
				if(key.equals("ref"))
				{
					refFn = val;
				}
				else if(key.equals("vcf"))
				{
					vcfFn = val;
				}
				else if(key.equals("output"))
				{
					outFn = val;
				}
			}
		}
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		// Check that reference file provided is valid
		if(!new File(refFn).exists())
		{
			usage();
			System.err.println("Reference file does not exist: " + refFn);
			System.exit(1);
		}
		
		// Check that variant file provided is valid
		if(!new File(vcfFn).exists())
		{
			usage();
			System.err.println("VCF does not exist: " + vcfFn);
			System.exit(1);
		}
		
		// Read the variants from the VCF file
		HashMap<String, TreeSet<Variant>> vars = new HashMap<String, TreeSet<Variant>>();
		Scanner vcfInput = new Scanner(new FileInputStream(new File(vcfFn)));
		while(vcfInput.hasNext())
		{
			String line = vcfInput.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			if(line.startsWith("#"))
			{
				continue;
			}
			String[] tokens = line.split("\t");
			String chr = tokens[0];
			int pos = Integer.parseInt(tokens[1]);
			String ref = tokens[3];
			String alt = tokens[4];
			
			if(tokens.length >= 10 && tokens[9].length() >= 3)
			{
				String genotype = tokens[9].substring(0, 3);
				if(genotype.equals("0|0") || genotype.equals("0/0"))
				{
					continue;
				}
				if(ignoreHeterozygous)
				{
					if(genotype.equals("1|0") || genotype.equals("1/0"))
					{
						continue;
					}
					if(genotype.equals("0|1") || genotype.equals("0/1"))
					{
						continue;
					}
				}
			}
			if(!vars.containsKey(chr))
			{
				vars.put(chr, new TreeSet<Variant>());
			}
			vars.get(chr).add(new Variant(pos, ref, alt));
		}
		vcfInput.close();
		
		// Read through the genome and output the updated version
		Scanner refInput = new Scanner(new FileInputStream(new File(refFn)));
		PrintWriter out = new PrintWriter(new File(outFn));
		
		String chrNameLine = "";
		StringBuilder contigSeq = new StringBuilder("");
		while(refInput.hasNext())
		{
			String line = refInput.nextLine();
			if(line.startsWith(">"))
			{
				if(chrNameLine.length() > 0)
				{
					String newContig = updateContig(contigSeq.toString(), chrNameLine, vars);
					printContig(newContig, chrNameLine, out);
				}
				chrNameLine = line;
				contigSeq = new StringBuilder("");
			}
			else
			{
				contigSeq.append(line);
			}
		}
		String newContig = updateContig(contigSeq.toString(), chrNameLine, vars);
		printContig(newContig, chrNameLine, out);
		
		refInput.close();
		out.close();
	}
	
	/*
	 * Prints a contig to the given output stream
	 */
	static void printContig(String s, String chrNameLine, PrintWriter out)
	{
		out.println(chrNameLine);
		int n = s.length();
		for(int i = 0; i<n; i++)
		{
			out.print(s.charAt(i));
			if(i == n-1 || i%charsPerLine == charsPerLine-1)
			{
				out.println();
			}
		}
	}
	
	/*
	 *  Integrates variants into a contig
	 */
	static String updateContig(String s, String chrNameLine, HashMap<String, TreeSet<Variant>> vars)
	{
		String name = chrNameLine.split(" ")[0].substring(1);
		
		TreeSet<Variant> varSet = vars.get(name);
		
		int n = s.length();
		int nextVariant = varSet.size() > 0 ? (varSet.first().pos - 1) : n;
		
		StringBuilder res = new StringBuilder("");
		
		for(int i = 0; i<n; i++)
		{
			while(nextVariant < i)
			{
				varSet.pollFirst();
				nextVariant = varSet.size() > 0 ? (varSet.first().pos - 1) : n;
			}
			if(i == nextVariant)
			{
				Variant var = varSet.pollFirst();
				if(!var.ref.equalsIgnoreCase(s.substring(i, i + var.ref.length())))
				{
					System.err.println("Ref sequence at position i (" + 
							s.substring(i, i + var.ref.length()) + ") doesn't match variant REF sequence: " + var.ref);
				}
				i += var.ref.length() - 1;
				res.append(var.alt);
				nextVariant = varSet.size() > 0 ? (varSet.first().pos - 1) : n;
			}
			else
			{
				res.append(s.charAt(i));
			}
		}
		return res.toString();
	}
	
	/*
	 * A variant to be integrated into the genome
	 */
	static class Variant implements Comparable<Variant>
	{
		int pos;
		String ref, alt;
		Variant(int pos, String ref, String alt)
		{
			this.pos = pos;
			this.ref = ref;
			this.alt = alt;
		}
		public int compareTo(Variant o)
		{
			return pos - o.pos;
		}
	}
}
