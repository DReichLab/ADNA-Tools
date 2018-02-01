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
	public static void main(String[] args) throws ParseException, IOException, java.text.ParseException {
		final String versionOption = "version";
		String command;
		String [] remainingArgs = new String[0];
		
		if(args.length == 0)
			command = versionOption;
		else {
			command = args[0];
			remainingArgs = Arrays.copyOfRange(args, 1, args.length);
		}
		switch(command.toLowerCase()){
		case "barcodecount":
			BarcodeCount.main(remainingArgs);
			break;
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
		case "assignreadgroups":
			AssignReadGroups.main(remainingArgs);
			break;
		case versionOption:
			System.out.println(versionString());
			break;
		default:
			throw new IllegalArgumentException("Unknown program: " + command);
		}
	}
	
	public static String versionString() {
		String version = Driver.class.getPackage().getImplementationVersion();
		if (version == null)
			version = "development";
		return version;
	}

}
