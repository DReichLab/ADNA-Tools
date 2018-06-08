package adnascreen;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.Stream;

/**
 * For each experiment key (2 indices and 2 barcodes), keep a SampleCounter of counts with labels
 *
 */
public class SampleSetsCounter {
	private long raw;
	private Map<IndexAndBarcodeKey, SampleCounter> sets;
	
	public SampleSetsCounter(){
		raw = 0;
		sets = new HashMap<IndexAndBarcodeKey, SampleCounter>();
	}
	
	public SampleSetsCounter(String s){
		this();

		String [] lines = s.split("\n");
		// first line is raw count
		raw = Long.valueOf(lines[0]);
		for(int n = 1; n < lines.length; n++){
			addKeyLine(lines[n]);
		}
	}
	
	public SampleSetsCounter(SampleSetsCounter other){
		this();
		this.raw = other.raw;
		other.sets.forEach((key, sample)->{
			this.sets.put(key, new SampleCounter(sample));
		});
	}
	
	public SampleSetsCounter(File input) throws IOException{
		this();
		try(BufferedReader reader = new BufferedReader(new FileReader(input))){
			raw = Long.valueOf(reader.readLine());
			reader.lines().forEach((String s) -> {
				addKeyLine(s);
			});
		}
	}
	
	/**
	 * Take a line from the String or file representation and add the 
	 * @param s tab separated key, and (label, count) pairs
	 */
	private void addKeyLine(String s){
		String [] fields = s.split("\t");
		// leftmost field is key
		IndexAndBarcodeKey key = new IndexAndBarcodeKey(fields[0]);
		// followed by pairs of labels with counts
		for(int i = 1; i + 1 < fields.length; i += 2){
			String label = fields[i];
			int count = Integer.valueOf(fields[i + 1]);
			this.add(key, label, count);
		}
	}
	
	public long increment(){
		return ++raw;
	}
	
	public long add(long value){
		raw += value;
		return raw;
	}
	
	public int increment(IndexAndBarcodeKey key, String label){
		SampleCounter sample = sets.get(key);
		if(sample == null){
			sample = new SampleCounter();
			sets.put(key, sample);
		}
		return sample.increment(label);
	}
	
	public int add(IndexAndBarcodeKey key, String label, int value){
		SampleCounter sample = sets.get(key);
		if(sample == null){
			sample = new SampleCounter();
			sets.put(key, sample);
		}
		return sample.add(label, value);
	}
	
	public SampleCounter get(IndexAndBarcodeKey key){
		return sets.get(key);
	}
	
	public int get(IndexAndBarcodeKey key, String label){
		SampleCounter s = sets.get(key);
		int value = (s == null ? 0 : s.get(label));
		return value;
	}
	
	public void combine(SampleSetsCounter other){
		this.raw += other.raw;
		other.sets.forEach((key, sample)->{
			// combine if exists
			if(sets.containsKey(key)){
				sets.get(key).combine(sample);
			} else { // append if it does not
				sets.put(key, new SampleCounter(sample));
			}
		});
	}
	
	public String toStringSorted(final String label){
		// sort by descending count in this label
		Stream<Map.Entry<IndexAndBarcodeKey, SampleCounter>> sortedByCount = sets.entrySet().stream().sorted(
				new Comparator<Map.Entry<IndexAndBarcodeKey, SampleCounter>>(){
					public int compare(Map.Entry<IndexAndBarcodeKey, SampleCounter> a, Map.Entry<IndexAndBarcodeKey, SampleCounter> b){
						Integer a_count = a.getValue().get(label);
						Integer b_count = b.getValue().get(label);
						return a_count.compareTo(b_count);
					}
				}.reversed());
		// convert this sorted stream into a string
		StringBuilder builder = new StringBuilder();
		builder.append(raw);
		builder.append('\n');
		sortedByCount.forEachOrdered((pair)->{
			IndexAndBarcodeKey key = pair.getKey();
			SampleCounter sample = pair.getValue();
			builder.append(key.toString());
			builder.append('\t');
			builder.append(sample.toString());
			builder.append('\n');
		});
		builder.deleteCharAt(builder.length()-1); // remove last newline
		return builder.toString();
	}
	
	@Override
	public String toString(){
		StringBuilder builder = new StringBuilder();
		builder.append(raw);
		builder.append('\n');
		sets.forEach((key, sample)->{
			builder.append(key.toString());
			builder.append('\t');
			builder.append(sample.toString());
			builder.append('\n');
		});
		builder.deleteCharAt(builder.length()-1); // remove last newline
		return builder.toString();
	}
	
	@Override
	public boolean equals(Object x){
		if(x instanceof SampleSetsCounter){
			SampleSetsCounter other = (SampleSetsCounter) x;
			// check in both directions
			for(IndexAndBarcodeKey key : sets.keySet()){
				if(!sets.get(key).equals(other.get(key))){
					return false;
				}
			}
			for(IndexAndBarcodeKey key : other.sets.keySet()){
				if(!other.get(key).equals(sets.get(key))){
					return false;
				}
			}
			return (this.raw == other.raw);
		}
		return false;
	}
}
