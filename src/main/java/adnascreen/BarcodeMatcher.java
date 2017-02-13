package adnascreen;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * This class is similar to the IndexMatcher, but for sets of barcodes. 
 * All barcodes within a set are considered equivalent. 
 * When searching, we return a string representing the entire set. 
 * @author mmah
 *
 */
public class BarcodeMatcher {
	private Set<String> referenceSets;
	private Map<DNASequence, String> cache;
	private int maxHammingDistance;
	private int barcodeLength = -1;
	
	public BarcodeMatcher(){
		referenceSets = new HashSet<String>();
		cache = new HashMap<DNASequence, String>();
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
	 * @param filename each line of file contains colon-delimited 
	 * @throws IOException
	 */
	public void loadFile(String filename) throws IOException {
		try(InputStream referenceStream = new FileInputStream(filename);
				BufferedReader reader = new BufferedReader(new InputStreamReader(referenceStream));
				){
			String line;
			while((line = reader.readLine()) != null){
				// each line is a barcode set
				String barcodeSetString = line.trim();
				referenceSets.add(barcodeSetString);
			}
			cache.clear();
			seedCache();
		}
	}
	
	public void addReferenceSet(String referenceSet){
		this.referenceSets.add(referenceSet);
		cache.clear();
		seedCache();
	}
	
	private void seedCache(){
		for(String barcodeSetString : referenceSets){
			String [] barcodeStrings = barcodeSetString.split(":");
			// check that barcodes are of the same length
			if(barcodeStrings.length > 0){
				for(String barcode : barcodeStrings){
					if(barcodeLength == -1) // set initial length once
						barcodeLength = barcode.length();
					if(barcode.length() != barcodeLength) // check all lengths against the first
						throw new IllegalArgumentException("barcode length mismatch");
					cache.put(new DNASequence(barcode), barcodeSetString);
				}
			}
		}
	}
	
	private String linearSearch(DNASequence query){
		int minDistance = maxHammingDistance + 1;
		String bestSet = null;
		for(String barcodeSetString : referenceSets){
			String [] barcodeStrings = barcodeSetString.split(":");
			for(String barcodeString : barcodeStrings){
				DNASequence barcode = new DNASequence(barcodeString);
				int distance = query.hammingDistance(barcode);
				if(distance < minDistance){
					minDistance = distance;
					bestSet = barcodeSetString;
					if(distance == 0){
						return bestSet;
					}
				}
			}
		}
		return bestSet;
	}
	
	public String find(DNASequence query){
		String found = null;
		// is value in cache?
		found = cache.get(query);
		if(found == null){ // not in cache
			// perform linear search
			found = linearSearch(query);
			// cache value
			cache.put(query, found);
		}
		return found;
	}
	
	public void setMaxHammingDistance(int maxHammingDistance){
		if(maxHammingDistance < 0){
			throw new IllegalArgumentException();
		}
		if(maxHammingDistance < this.maxHammingDistance)
			cache.clear();
		this.maxHammingDistance = maxHammingDistance;
	}
}
