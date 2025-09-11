package adnatools;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;

public class FASTQFilter {

	public static void main(String[] args) throws IOException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("n", "length", true, "minimum length to retain");
		options.addOption("i", "input", true, "Input FASTQ.gz file. Otherwise read from standard input");
		options.addOption("o", "output", true, "Output FASTQ.gz file. Otherwise write to standard output");
		
		try {
			CommandLine commandLine	= parser.parse(options, args);
			int minimum_length = Integer.valueOf(commandLine.getOptionValue("length"));
			// construct input stream
			FastqReader in;
			if (commandLine.hasOption("input")){
				String inputFilename = commandLine.getOptionValue("input");
				in = new FastqReader(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inputFilename)))));
			} else {
				in = new FastqReader(new BufferedReader(new InputStreamReader(System.in)));
			}
			// construct output stream
			PrintWriter out;
			if (commandLine.hasOption("output")){
				String outputFilename = commandLine.getOptionValue("output");
				out = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outputFilename)));
			} else {
				out = new PrintWriter(System.out);
			}
			
			filter(in, out, minimum_length);
			in.close();
			out.close();
		} catch (ParseException e) {
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp("FASTQFilter", options);
			System.err.println(e.getLocalizedMessage());
		}
	}
	
	public static void filter(FastqReader in, PrintWriter out, int minimum_length) {
		while(in.hasNext()) {
			FastqRecord record = in.next();
			if(record.getReadLength() >= minimum_length) {
				out.write(record.toFastQString());
				out.write('\n');
			}
		}
	}

}
