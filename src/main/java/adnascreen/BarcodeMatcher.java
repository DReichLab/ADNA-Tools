package adnascreen;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Optional;
import java.util.Set;

/**
 * Find a best match for a query within a specified tolerance (maximum Hamming distance)
 * within a known reference set. 
 * Possible future improvements:
 * - comparison based on probability rather than Hamming distance
 * @author mmah
 *
 */
public class BarcodeMatcher {
	private Map<DNASequence, String> referenceBarcodeToLabel;
	private Map<DNASequence, Optional<String> >cache;
	private Set<String> labels;
	private int maxHammingDistance;
	private int barcodeLength = -1;
	public final static char INDEX_DELIMITER = '.';
	
	public BarcodeMatcher(){
		referenceBarcodeToLabel = new HashMap<DNASequence, String>();
		cache = new HashMap<DNASequence, Optional<String> >();
		labels = new HashSet<String>();
	}
	
	public BarcodeMatcher(String filename, int maxHammingDistance) throws IOException{
		this();
		this.maxHammingDistance = maxHammingDistance;
		loadFile(filename);
	}
	
	public int getBarcodeLength(){
		return barcodeLength;
	}

	/**
	 * 
	 * @param filename each line of file contains colon-delimited barcodes followed by whitespace and label
	 * @throws IOException
	 */
	public void loadFile(String filename) throws IOException {
		try(InputStream referenceStream = new FileInputStream(filename);
				BufferedReader reader = new BufferedReader(new InputStreamReader(referenceStream));
				){
			String line;
			while((line = reader.readLine()) != null){
				// each line is a barcode set string, with whitespace and label
				String [] elements = line.trim().split("\\s+");
				String barcodeSetString = elements[0];
				String barcodeLabel = elements[1];
				addReferenceSet(barcodeSetString, barcodeLabel, false);
			}
			cache.clear();
			seedCache();
		}
	}
	
	/**
	 * 
	 * @param barcodeSetString colon-delimited set of DNA barcode strings
	 * @param label unique descriptor for set
	 * @param clearCaches When adding multiple references sets, this should be set to false, 
	 * then cache should be called manually. 
	 */
	private void addReferenceSet(String barcodeSetString, String label, boolean clearCaches){
		String [] barcodeStrings = barcodeSetString.split(":");
		// check that barcodes are of the same length
		if(barcodeStrings.length > 0){
			// enforce unique labels
			if(labels.contains(label))
				throw new IllegalArgumentException("labels must be unique");
			if(label.contains(String.valueOf(INDEX_DELIMITER)))
				throw new IllegalArgumentException("Barcode labels can not contain delimiter " + INDEX_DELIMITER);
			labels.add(label);
			// add reference sets for this label
			int index = 1;
			for(String barcode : barcodeStrings){
				if(barcodeLength == -1) // set initial length once
					barcodeLength = barcode.length();
				if(barcode.length() != barcodeLength) // check all lengths against the first
					throw new IllegalArgumentException("barcode length mismatch");
				// if there is more than one barcode in this set, add index within set 
				String augmentedLabel = (barcodeStrings.length > 1) ? label + INDEX_DELIMITER + index++ : label;
				DNASequence barcodeSequence = new DNASequence(barcode);
				if(referenceBarcodeToLabel.containsKey(barcodeSequence))
					throw new IllegalArgumentException("barcodes must be unique");
				else
					referenceBarcodeToLabel.put(barcodeSequence, augmentedLabel);
			}
		}
		if(clearCaches){
			cache.clear();
			seedCache();
		}
	}
	
	/**
	 * 
	 * @param barcodeSetString colon-delimited set of DNA barcode strings
	 * @param label unique descriptor for set
	 */
	public void addReferenceSet(String barcodeSetString, String label){
		addReferenceSet(barcodeSetString, label, true);
	}
	
	private void seedCache(){
		for(DNASequence barcode : referenceBarcodeToLabel.keySet()){
			String label = referenceBarcodeToLabel.get(barcode);
			cache.put(barcode, Optional.of(label));
		}
	}
	
	/**
	 * 
	 * @param query
	 * @return label matching query within maxHammingDistance, or null if no match
	 */
	private String linearSearch(DNASequence query){
		int bestDistance = maxHammingDistance + 1;
		String bestLabel = null;
		for(DNASequence barcode : referenceBarcodeToLabel.keySet()){
			int distance = query.hammingDistance(barcode);
			if(distance < bestDistance){
				bestDistance = distance;
				bestLabel = referenceBarcodeToLabel.get(barcode);
				if(distance == 0){
					return bestLabel;
				}
			}
		}
		return bestLabel;
	}
	
	/**
	 * 
	 * @param query
	 * @return label matching or nearly matching query
	 */
	public String find(DNASequence query){
		// use Optional to differentiate between 
		// 1. query not present in cache
		// 2. cache knows there is no value for this query
		Optional<String> found = null;
		// is value in cache?
		found = cache.get(query);
		if(found == null){ // not in cache
			// perform linear search
			String label = linearSearch(query);
			found = Optional.ofNullable(label) ;
			// cache possibly null value
			cache.put(query, found);
		}
		return found.orElse(null);
	}
	
	/**
	 * This invalidates any cache, so it should be called before inserting reference sets. 
	 * @param maxHammingDistance maximum Hamming distance allowed for matches
	 */
	public void setMaxHammingDistance(int maxHammingDistance){
		if(maxHammingDistance < 0){
			throw new IllegalArgumentException();
		}
		if(maxHammingDistance < this.maxHammingDistance)
			cache.clear();
		this.maxHammingDistance = maxHammingDistance;
	}
}
