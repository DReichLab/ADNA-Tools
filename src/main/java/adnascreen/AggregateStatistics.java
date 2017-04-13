package adnascreen;

import java.io.File;
import java.io.IOException;

public class AggregateStatistics {
	/**
	 * Take a series of SampleCounter outputs and aggregate them
	 * @param args
	 */
	public static void main(String [] args) throws IOException {
		// read statistics from set of files
		SampleSetsCounter statistics = new SampleSetsCounter();
		for(String filename : args){
			File f = new File(filename);
			SampleSetsCounter current = new SampleSetsCounter(f);
			statistics.combine(current);
		}
		// output combined statistics
		System.out.println(statistics.toStringSorted("raw"));
	}
}
