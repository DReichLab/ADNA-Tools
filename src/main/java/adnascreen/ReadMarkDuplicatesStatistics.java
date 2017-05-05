package adnascreen;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

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
	
	public static IndexAndBarcodeKey keyFromFilename(String filenameFullPath){
		String [] path = filenameFullPath.split("/");
		String filename = path[path.length - 1]; // filename is at the end of the path
		String [] parts = filename.split("\\.");
		String filenameOnly = parts[0]; // key is at beginning of filename
		IndexAndBarcodeKey key = new IndexAndBarcodeKey(filenameOnly); 
		return key;
	}

	public static void main(String args[]) throws IOException{
		String filenameFullPath = args[0];
		IndexAndBarcodeKey key = keyFromFilename(filenameFullPath);
		// read mark duplicates output
		SampleSetsCounter stats = new SampleSetsCounter();
		int numDuplicates = readMarkDuplicateStatistics(filenameFullPath);
		stats.add(key, DUPLICATES, numDuplicates);
		System.out.println(stats.toString());
	}
}
