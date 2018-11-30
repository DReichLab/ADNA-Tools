package adnascreen;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math3.stat.Frequency;
import org.json.JSONArray;
import org.json.JSONObject;

import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

/**
 * Count the number of alignments for specified chromosome targets
 *
 */
public class SAMStats {
	private Frequency lengthHistogram;
	// each alignment may map to multiple targets
	private Map<String, Set<String> > referenceNamesToTargetClasses;
	SampleSetsCounter counter;

	public SAMStats(String keyString, String filename, JSONObject targets, int minimumMappingQuality) throws IOException{
		lengthHistogram = new Frequency();
		counter = new SampleSetsCounter();
		referenceNamesToTargetClasses = new HashMap<String, Set<String> >();
		
		// construct data structures for counting
		for(String targetClass : targets.keySet()){
			JSONArray array = targets.optJSONArray(targetClass);			
			if(array != null){
				for(int i = 0; i < array.length(); i++){
					String referenceName = array.getString(i);
					addReferenceToTargetClass(referenceName, targetClass);
				}
			} else {
				String referenceName = targets.getString(targetClass);
				addReferenceToTargetClass(referenceName, targetClass);
			} 
		}
		
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
				){
			//SAMFileHeader header = reader.getFileHeader();
			//SAMSequenceDictionary currentAlignmentReference = header.getSequenceDictionary();

			SAMRecordIterator i = reader.iterator();
			while(i.hasNext()){
				// iterate through alignments
				try{
					SAMRecord record = i.next();
					if(!record.getReadUnmappedFlag()){ // read is mapped
						int mappingQuality = record.getMappingQuality();
						if(mappingQuality >= minimumMappingQuality){
							// record length
							int readLength = record.getReadLength();
							lengthHistogram.addValue(readLength);
							// count alignments to reference sequences 
							String referenceName = record.getReferenceName();
							Set<String> targetClassesToCount = referenceNamesToTargetClasses.get(referenceName);
							if(targetClassesToCount != null){
								for(String targetClass : targetClassesToCount){
									counter.increment(keyString, targetClass);
									String targetCoverageLabel = targetClass + "-coverageLength";
									counter.add(keyString, targetCoverageLabel, readLength);
								}
							}
							//System.out.println(referenceName + "\t" + targetClass);
						}
					}
				} catch (SAMFormatException e){
					System.err.println(e);
					// ignore this record and continue to the next
				}
			}
		}
	}
	
	private void addReferenceToTargetClass(String reference, String targetClass){
		if(!referenceNamesToTargetClasses.containsKey(reference)){
			referenceNamesToTargetClasses.put(reference, new HashSet<String>() );
		}
		Set<String> matchingTargets = referenceNamesToTargetClasses.get(reference);
		matchingTargets.add(targetClass);
	}
	
	public String toString(){
		return counter.toString();
	}
	
	public static void main(String []args) throws IOException, ParseException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("f", "filename", true, "input SAM/BAM filename");
		options.addRequiredOption("t", "targets", true, "JSON object for targets: key is label, value is reference name or array of reference names");
		options.addOption("l", "length-histogram", true, "filename for length distribution histogram");
		options.addOption("q", "mapping-quality", true, "minimum Phred score mapping quality");
		CommandLine commandLine	= parser.parse( options, args );
		
		String filename = commandLine.getOptionValue('f');
		JSONObject targets = new JSONObject(commandLine.getOptionValue('t'));
		String histogramFilename = commandLine.getOptionValue('l');
		int minimumMappingQuality = Integer.valueOf(commandLine.getOptionValue('q', "0"));
		
		String keyString = ReadMarkDuplicatesStatistics.keyFromFilename(filename);
		SAMStats stats = new SAMStats(keyString, filename, targets, minimumMappingQuality);
		System.out.println(stats.toString());
		
		if(histogramFilename != null){
			try(
					FileWriter histogramFile = new FileWriter(histogramFilename);
					PrintWriter w = new PrintWriter(histogramFile)
					){
				w.println(stats.lengthHistogram.toString());
			}
		}
	}
}
