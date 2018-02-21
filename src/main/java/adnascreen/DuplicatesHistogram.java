package adnascreen;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.lang3.mutable.MutableInt;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class DuplicatesHistogram {
	int maximumDepth = 0;
	HashMap<Integer, MutableInt> occurrencesByDepth;
	
	DuplicatesHistogram(){
		occurrencesByDepth = new HashMap<Integer, MutableInt>();
	}
	
	DuplicatesHistogram(String inputFilename) throws IOException{
		occurrencesByDepth = new HashMap<Integer, MutableInt>();
		uniqueReadHistogram(inputFilename);
	}
	
	public static void main(String [] args) throws ParseException, IOException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input_BAM", true, "Input BAM filename");
		CommandLine commandLine	= parser.parse(options, args);
		
		String inputFilename = commandLine.getOptionValue('i');
		DuplicatesHistogram duplicatesHistogram = new DuplicatesHistogram(inputFilename);
		int[] histogram = duplicatesHistogram.getHistogram();
		for(int i = 0; i < histogram.length; i++) {
			System.out.println( (i+1) + "\t" + histogram[i]);
		}
	}
	
	private void uniqueReadHistogram(String inputFilename) throws IOException {		
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename)));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
				){
			SAMRecordIterator i = reader.iterator();
			HashMap<String, MutableInt> countsAtThisPosition = new HashMap<String, MutableInt>();
			String currentReference = null;
			int currentAlignmentStart = 0;
			while(i.hasNext()){
				// iterate through alignments
				try{
					SAMRecord record = i.next();
					 
					if(!record.getReadUnmappedFlag()) {
						String reference = record.getReferenceName();
						int start = record.getAlignmentStart();
						// if we are at a new start position, record contents of the current position and reset 
						if(!reference.equals(currentReference) || currentAlignmentStart != start) {
							updateHistogramWithCurrentStartPosition(countsAtThisPosition);
							countsAtThisPosition.clear();
						}
						// at this start position, count
						currentReference = reference;
						currentAlignmentStart = start;
						String duplicateKey = reference + "_" + start + "_" + record.getStringAttribute("XD");
						if(!countsAtThisPosition.containsKey(duplicateKey)) {
							countsAtThisPosition.put(duplicateKey, new MutableInt(0));
						}
						MutableInt counter = countsAtThisPosition.get(duplicateKey);
						counter.increment();
					}
				}
				catch(Exception e){
					System.err.println(e.toString());
				}
			}
			updateHistogramWithCurrentStartPosition(countsAtThisPosition);
			countsAtThisPosition.clear();
		}
	}
	
	private void updateHistogramWithCurrentStartPosition(HashMap<String, MutableInt> countsAtThisPosition) {
		for(String key : countsAtThisPosition.keySet()) {
			int depth = countsAtThisPosition.get(key).getValue();
			maximumDepth = Math.max(depth, maximumDepth);
			if(!occurrencesByDepth.containsKey(depth)) {
				occurrencesByDepth.put(depth, new MutableInt(0));
			}
			MutableInt count = occurrencesByDepth.get(depth);
			count.increment();
		}
	}
	
	public int[] getHistogram() {
		// populate histogram
		int [] histogram = new int[maximumDepth];
		for(int depth = 1; depth <= maximumDepth; depth++) {
			histogram[depth-1] = occurrencesByDepth.getOrDefault(depth, new MutableInt(0)).getValue();
		}
		return histogram;
	}
}
