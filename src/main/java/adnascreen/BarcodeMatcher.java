package adnascreen;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.SortedSet;
import java.util.TreeSet;

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
	private Map<String, Integer> labelToBarcodeLength;
	private SortedSet<Integer> barcodeLengthsSet;
	private List<Integer> barcodeLengths;
	private int maxHammingDistance;
	public final static char INDEX_DELIMITER = '.';
	public final static char BARCODE_DELIMITER = ':';
	private Map<String, DNASequence> labelToBarcode;
	
	public BarcodeMatcher(){
		referenceBarcodeToLabel = new HashMap<DNASequence, String>();
		cache = new HashMap<DNASequence, Optional<String> >();
		labelToBarcodeLength = new HashMap<String, Integer>();
		barcodeLengthsSet = new TreeSet<Integer>();
		labelToBarcode = new HashMap<String, DNASequence>();
	}
	
	public BarcodeMatcher(String filename, int maxHammingDistance) throws IOException{
		this();
		this.maxHammingDistance = maxHammingDistance;
		loadFile(filename);
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
			barcodeLengths = null;
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
		String [] barcodeStrings = barcodeSetString.toUpperCase().split(String.valueOf(BARCODE_DELIMITER));
		// check that barcodes are of the same length
		if(barcodeStrings.length > 0){
			int barcodeLength = -1;
			// enforce unique labels
			if(labelToBarcodeLength.containsKey(label))
				throw new IllegalArgumentException("labels must be unique");
			if(label.contains(String.valueOf(INDEX_DELIMITER)))
				throw new IllegalArgumentException("Barcode labels can not contain delimiter " + INDEX_DELIMITER);
			// add reference sets for this label
			int index = 1;
			for(String barcode : barcodeStrings){
				if(barcodeLength == -1)
					barcodeLength = barcode.length();
				if(barcode.length() != barcodeLength) // check all lengths against the first
					throw new IllegalArgumentException("barcode length mismatch");
				// Fill out barcode to label map
				// if there is more than one barcode in this set, add index within set 
				String augmentedLabel = (barcodeStrings.length > 1) ? label + INDEX_DELIMITER + index++ : label;
				DNASequence barcodeSequence = new DNASequence(barcode);
				if(referenceBarcodeToLabel.containsKey(barcodeSequence))
					throw new IllegalArgumentException("barcodes must be unique");
				else
					referenceBarcodeToLabel.put(barcodeSequence, augmentedLabel);
				// Fill out label to barcode map (reverse of previous)
				labelToBarcode.put(augmentedLabel, barcodeSequence);
			}
			if(barcodeLength > 0){
				labelToBarcodeLength.put(label, barcodeLength);
				barcodeLengthsSet.add(barcodeLength);
			}
		}
		if(clearCaches){
			cache.clear();
			seedCache();
			barcodeLengths = null;
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
			if(query.length() == barcode.length()){
				int distance = query.hammingDistance(barcode);
				if(distance < bestDistance){
					bestDistance = distance;
					bestLabel = referenceBarcodeToLabel.get(barcode);
					if(distance == 0){
						return bestLabel;
					}
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
	
	/**
	 * 
	 * @param label
	 * @return length of barcode(s) corresponding to this label, 0 barcode label is not found
	 */
	public int getBarcodeLength(String label){
		String flattenedLabel = IndexAndBarcodeKey.flatten(label);
		Integer length = labelToBarcodeLength.get(flattenedLabel);
		return (length != null) ? length : 0;
	}
	
	/**
	 * This is a convenience function to parse the barcode portion of 
	 * an IndexAndBarcodeKey to determine barcode length
	 * @param pairLabel
	 * @return barcode length, or 0 if barcode label is not found
	 */
	public int getBarcodePairLength(String pairLabel){
		String[] maxBarcodeLabels = pairLabel.split(String.valueOf(IndexAndBarcodeKey.FIELD_SEPARATOR));
		int barcode1Length = getBarcodeLength(maxBarcodeLabels[0]);
		if(maxBarcodeLabels.length > 1){
			int barcode2Length = getBarcodeLength(maxBarcodeLabels[1]);
			if(barcode1Length != barcode2Length)
				throw new IllegalArgumentException("barcode lengths do not match");
		}
		return barcode1Length;
	}
	
	private void cacheBarcodeLengths(){
		List<Integer> lengths = new ArrayList<Integer>(barcodeLengthsSet);
		Collections.reverse(lengths);
		barcodeLengths = lengths;
	}
	
	/**
	 * Get a list of barcode lengths sorted in descending order.
	 * Do NOT modify this list without copying it. 
	 * @return 
	 */
	public List<Integer> getBarcodeLengths(){
		if(barcodeLengths == null)
			cacheBarcodeLengths();
		return barcodeLengths;
	}
	
	public DNASequence getBarcode(String label) {
		return labelToBarcode.get(label);
	}
}
