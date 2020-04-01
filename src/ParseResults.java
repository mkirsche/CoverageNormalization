import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.Scanner;
import java.util.TreeMap;
public class ParseResults {
	static boolean homo_only = false;
	static String OUTDIR = "output";
public static void main(String[] args) throws Exception
{
	int cov = 200;
	if(args.length > 0)
	{
		cov = Integer.parseInt(args[0]);
		OUTDIR = args[1];
		if(args.length > 2)
		{
			if(args[2].equalsIgnoreCase("--homozygous"))
			{
				homo_only = true;
			}
		}
	}
	String fn = OUTDIR+"/vcfs_"+cov+".txt";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	TreeMap<Integer, Integer> fs = new TreeMap<Integer, Integer>();
	int count = 0;
	while(input.hasNext())
	{
		String line = input.nextLine();
		String[] tokens = line.split("\t");
		if(tokens.length >= 10)
		{
			if(!homo_only || (tokens[9].startsWith("1/1") || tokens[9].startsWith("1|1")))
			{
				int position = Integer.parseInt(tokens[1]);
				fs.put(position, 1 + fs.getOrDefault(position, 0));
			}
		}
		else if(tokens.length == 1) count++;
	}
	String ofn = OUTDIR+"/counts_"+cov+".txt";
	if(homo_only) ofn = OUTDIR+"/counts_homo_"+cov+".txt";
	PrintWriter out = new PrintWriter(new File(ofn));
	for(int x : fs.keySet())
	{
		out.println(x+" "+(fs.get(x)*1.0/count));
	}
	out.close();
	input.close();
}
}
