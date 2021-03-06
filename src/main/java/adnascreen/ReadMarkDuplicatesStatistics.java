package adnascreen;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class ReadMarkDuplicatesStatistics {
	public static final String DUPLICATES = "duplicates";
	/**
	 * Read the number of unpaired duplicates from Picard output
	 * @param filename for Picard MarkDuplicates output
	 * @return number of duplicates of unpaired (merged) reads
	 * @throws FileNotFoundException
	 * @throws IOException
	 */
	public static int readMarkDuplicateStatistics(String filename) throws FileNotFoundException, IOException{
		try(BufferedReader f = new BufferedReader(new FileReader(filename))){
			String line;
			// ignore comment lines beginning with # and empty lines
			do{
				line = f.readLine();
			} while(line != null && (line.startsWith("#") || line.length() == 0));
			if(line == null) // no line indicates there are no duplicates
				return 0;
			String headerLine = line;
			String statsLine = f.readLine();

			String[] headers = headerLine.split("\t");
			String[] stats = statsLine.split("\t");
			for(int n = 0; n < headers.length; n++){
				if(headers[n].equals("UNPAIRED_READ_DUPLICATES") && n < stats.length)
					return Integer.valueOf(stats[n]);
			}
		}
		return -1;	
	}
	
	public static String keyFromFilename(String filenameFullPath){
		String [] path = filenameFullPath.split("/");
		String filename = path[path.length - 1]; // filename is at the end of the path
		
		// remove extenions after .bam
		String filenameOnly;
		int bamPosition = filename.indexOf(".bam");
		if(bamPosition > 0) {
			filenameOnly = filename.substring(0, bamPosition);
		}
		else
			filenameOnly = filename;
		return filenameOnly;
	}

	public static void main(String args[]) throws IOException, ParseException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addOption("l", "label", true, "label for duplicates in statistics");
		CommandLine commandLine	= parser.parse(options, args);
		
		String label = commandLine.getOptionValue('l', DUPLICATES);
		
		String filenameFullPath = commandLine.getArgList().get(0);
		String keyString = keyFromFilename(filenameFullPath);
		// read mark duplicates output
		SampleSetsCounter stats = new SampleSetsCounter();
		int numDuplicates = readMarkDuplicateStatistics(filenameFullPath);
		stats.add(keyString, label, numDuplicates);
		System.out.println(stats.toString());
	}
}
