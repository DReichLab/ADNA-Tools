package adnascreen;

import java.io.BufferedInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

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
 * Count the number of alignments for specified targets
 *
 */
public class SAMStats {
	private Frequency lengthHistogram;
	private Map<String, String> referenceNamesToTargetClass;
	SampleSetsCounter counter;

	public SAMStats(String filename, JSONObject targets) throws IOException{
		lengthHistogram = new Frequency();
		IndexAndBarcodeKey key = ReadMarkDuplicatesStatistics.keyFromFilename(filename);
		counter = new SampleSetsCounter();
		referenceNamesToTargetClass = new HashMap<String, String>();
		
		for(String targetClass : targets.keySet()){
			JSONArray array = targets.optJSONArray(targetClass);
			if(array != null){
				for(int i = 0; i < array.length(); i++){
					String referenceName = array.getString(i);
					referenceNamesToTargetClass.put(referenceName, targetClass);
				}
			} else {
				String referenceName = targets.getString(targetClass);
				referenceNamesToTargetClass.put(referenceName, targetClass);
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
						// record length
						int readLength = record.getReadLength();
						lengthHistogram.addValue(readLength);
						// count alignments to reference sequences 
						String referenceName = record.getReferenceName();
						String targetClass = referenceNamesToTargetClass.get(referenceName);
						if(targetClass != null){
							counter.increment(key, targetClass);
						}
						//System.out.println(referenceName + "\t" + targetClass);
					}
				} catch (SAMFormatException e){
					System.err.println(e);
					// ignore this record and continue to the next
				}
			}
			System.out.println(counter);
		}
	}
	
	public static void main(String []args) throws IOException{
		String filename = args[0];
		JSONObject targets = new JSONObject(args[1]);
		String histogramFilename = null;
		if(args.length >= 2)
			histogramFilename = args[2];
		
		// parse the second argument as JSON
		// each key is a label
		// value is the reference name to associate with that label
		// or value is array of reference names to associate with label
		SAMStats stats = new SAMStats(filename, targets);
		
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
