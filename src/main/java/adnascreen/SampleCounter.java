package adnascreen;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.mutable.MutableInt;

public class SampleCounter {
	private Map<String, MutableInt> counters;
	private List<String> orderedLabels; // labels in order added to print in a specific order
	
	public SampleCounter(){
		counters = new HashMap<String, MutableInt>();
		orderedLabels = new LinkedList<String>();
	}
	
	// deep copy of other
	public SampleCounter(SampleCounter other){
		this();
		for(String label : other.orderedLabels){
			MutableInt value = new MutableInt(other.get(label));
			add(label, value.intValue());
		}
	}
	
	public SampleCounter(String serialized){
		this();
		String [] labelsAndValues = serialized.split("\t");
		for(int n = 0; n < labelsAndValues.length - 1; n += 2){
			String label = labelsAndValues[n];
			int value = Integer.valueOf(labelsAndValues[n+1]);
			add(label, value);
		}
	}
	
	public int increment(String label){
		MutableInt count = counters.get(label);
		if(count == null){
			count = new MutableInt(0);
			counters.put(label, count);
			orderedLabels.add(label);
		}
		return count.incrementAndGet();
	}
	
	public int add(String label, int value){
		MutableInt count = counters.get(label);
		if(count == null){
			count = new MutableInt(0);
			counters.put(label, count);
			orderedLabels.add(label);
		}
		return count.addAndGet(value);
	}
	
	public int get(String label){
		MutableInt count = counters.get(label);
		if(count != null){
			return count.intValue();
		}
		return 0;
	}
	
	/**
	 * Take counts from other SampleCounter and add them with counts from this one.
	 * Overwrites counts for this object in place. 
	 * @param other
	 */
	public void combine(SampleCounter other){
		for(String label : other.orderedLabels){
			int addend = other.get(label);
			add(label, addend);
		}
	}
	
	@Override
	public String toString(){
		StringBuilder builder = new StringBuilder();
		for(String label : orderedLabels){
				builder.append(label);	
				builder.append('\t');
				builder.append(counters.get(label));
				builder.append('\t');
		}
		builder.deleteCharAt(builder.length() - 1); // remove trailing tab
		return builder.toString();
	}
	
	@Override
	public boolean equals(Object x){
		if(x instanceof SampleCounter){
			SampleCounter other = (SampleCounter) x;
			// check for set inclusion both ways
			for(String label : other.orderedLabels){
				MutableInt otherValue = other.counters.get(label);
				MutableInt thisValue = counters.get(label);
				if(!otherValue.equals(thisValue)){
					return false;
				}
			}
			for(String label : orderedLabels){
				MutableInt otherValue = other.counters.get(label);
				MutableInt thisValue = counters.get(label);
				if(!thisValue.equals(otherValue)){
					return false;
				}
			}
			return true;
		}
		return false;
	}
}
