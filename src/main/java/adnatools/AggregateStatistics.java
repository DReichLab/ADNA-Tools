package adnatools;

import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.help.HelpFormatter;

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

		try {
			CommandLine commandLine	= parser.parse(options, args);

			String sortField = commandLine.getOptionValue('s', IndexAndBarcodeScreener.RAW);
			String NULL = null;
			String keyString = commandLine.getOptionValue('k', NULL);
			String label = commandLine.getOptionValue('l', NULL);

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
		} catch(ParseException e) {
			HelpFormatter formatter = HelpFormatter.builder().get();
			String commandDescription = "Restrict to only those reads meeting a post-mortem damage threshold.";
			formatter.printHelp(MethodHandles.lookup().lookupClass().getSimpleName(), commandDescription, options, Driver.submittingHelp, true);
			System.err.println(e.getLocalizedMessage());
			throw e;
		}
	}
}
