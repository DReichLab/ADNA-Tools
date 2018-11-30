package adnascreen;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class AggregateStatistics {
	/**
	 * Take a series of SampleCounter outputs and aggregate them
	 * @param args
	 */
	public static void main(String [] args) throws IOException, ParseException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addOption("s", "sort", true, "field to sort descending");
		options.addOption("k", "key", true, "sample key ");
		options.addOption("l", "label", true, "label to retrieve for key");
		CommandLine commandLine	= parser.parse(options, args);
		
		String sortField = commandLine.getOptionValue('s', IndexAndBarcodeScreener.RAW);
		String keyString = commandLine.getOptionValue('k', null);
		String label = commandLine.getOptionValue('l', null);
		
		// read statistics from set of files
		SampleSetsCounter statistics = new SampleSetsCounter();
		for(String filename : commandLine.getArgList()){
			File f = new File(filename);
			SampleSetsCounter current = new SampleSetsCounter(f);
			statistics.combine(current);
		}
		// retrieve one value if both key and label are specified
		if(keyString != null && label != null) {
			System.out.println(statistics.get(keyString, label));
		}
		else { // output combined statistics
			System.out.println(statistics.toStringSorted(sortField));
		}
	}
}
