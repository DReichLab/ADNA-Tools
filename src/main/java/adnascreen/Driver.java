package adnascreen;

import java.io.IOException;
import java.util.Arrays;

import org.apache.commons.cli.ParseException;

/**
 * This is a entry point to the ancient DNA screening programs. 
 * The first argument is used to determine which program to run. 
 * This program is passed the remaining arguments. 
 *
 */
public class Driver {

	public static void main(String[] args) throws ParseException, IOException {
		String command = args[0];
		String [] remainingArgs = Arrays.copyOfRange(args, 1, args.length);
		switch(command.toLowerCase()){
		case "indexandbarcodescreener":
			IndexAndBarcodeScreener.main(remainingArgs);
			break;
		case "aggregatestatistics":
			AggregateStatistics.main(remainingArgs);
			break;
		case "demultiplexsam":
			DemultiplexSAM.main(remainingArgs);
			break;
		case "readmarkduplicatesstatistics":
			ReadMarkDuplicatesStatistics.main(remainingArgs);
			break;
		case "samstats":
			SAMStats.main(remainingArgs);
			break;
		case "softclip":
			SoftClip.main(remainingArgs);
			break;
		default:
			throw new IllegalArgumentException("Unknown program: " + command);
		}
	}

}