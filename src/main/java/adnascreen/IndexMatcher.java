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
 * Find a best match for a query within a specified tolerance (maximum Hamming distance)
 * within a known reference set.  
 * Possible future improvements:
 * - comparison based on probability rather than Hamming distance
 * @author mmah
 *
 */
public class IndexMatcher {
	private Set<DNASequence> referenceSet;
	private Map<DNASequence, DNASequence> cache;
	private int maxHammingDistance = 0;
	
	public IndexMatcher(){
		referenceSet = new HashSet<DNASequence>();
		cache = new HashMap<DNASequence, DNASequence>();
	}
	
	/**
	 * 
	 * @param filename contains reference values, one per line
	 * @param maxHammingDistance
	 */
	public IndexMatcher(String filename, int maxHammingDistance) throws IOException{
		this();
		setMaxHammingDistance(maxHammingDistance);
		loadFile(filename);
	}
	
	public void loadFile(String filename) throws IOException{
		cache.clear();
		try(InputStream referenceStream = new FileInputStream(filename);
				BufferedReader reader = new BufferedReader(new InputStreamReader(referenceStream));
				){
			String line;
			int expectedLength = 0;
			while((line = reader.readLine()) != null){
				DNASequence referenceValue = new DNASequence(line.trim());
				int length = referenceValue.length();
				if(length > 0){
					if(expectedLength == 0){ // startup case
						expectedLength = length;
					}
					if(length == expectedLength){
						referenceSet.add(referenceValue);
					}
					else{
						throw new IllegalArgumentException("indices are different lengths");
					}
				}
				
			}
		}
	}
	
	/**
	 * Add a new reference to the set of possible 
	 * It is best to add all references before making queries because 
	 * adding a new reference invalidates the cache.
	 * @param reference
	 */
	public void addReference(DNASequence reference){
		this.referenceSet.add(reference);
		cache.clear();
	}
	
	private DNASequence linearSearch(DNASequence query){
		int minDistance = maxHammingDistance + 1;
		DNASequence best = null;
		for(DNASequence candidate : referenceSet){
			int distance = query.hammingDistance(candidate);
			if(distance < minDistance){
				minDistance = distance;
				best = candidate;
				if(distance == 0){
					return best;
				}
			}
		}
		return best;
	}
	
	public DNASequence find(DNASequence query){
		DNASequence found = null;
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
