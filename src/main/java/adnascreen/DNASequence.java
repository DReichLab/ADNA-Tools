package adnascreen;

import java.util.HashMap;
import java.util.Map;

public class DNASequence {
	private String sequence;
	
	public DNASequence(String s){
		sequence = s;
		// validation of string input for allowable characters
		for(int n = 0; n < s.length(); n++){
			char c = s.charAt(n);
			if(!reverseComplementMap.containsKey(c)){
				throw new IllegalArgumentException("Invalid DNASequence character " + c);
			}
		}
	}
	
	@Override
	public int hashCode(){
		return sequence.hashCode();
	}
	
	@Override
	public boolean equals(Object x){
		if(x instanceof DNASequence){
			DNASequence s = (DNASequence) x;
			return sequence.equals(s.sequence);
		}
		else
			return false;
	}
	
	public String toString(){
		return sequence;
	}
	
	public int length(){
		return sequence.length();
	}
	
	public char charAt(int index){
		return sequence.charAt(index);
	}

	public DNASequence reverseComplement(){
		StringBuilder b = new StringBuilder();
		for (int n = this.length() -1; n >= 0; n--){
			b.append(reverseComplement(sequence.charAt(n)));
		}
		return new DNASequence(b.toString());
	}
	
	public int hammingDistance(DNASequence s){
		int distance = 0;
		if(this.length() != s.length())
			throw new IllegalArgumentException();
		for(int n = 0; n < this.length(); n++){
			char c1 = charAt(n);
			char c2 = s.charAt(n);
			if(c1 != c2)
				distance += 1;
		}
		return distance;
	}
	
	protected static Map<Character, Character> reverseComplementMap;
	static{
		Map<Character, Character> map = new HashMap<Character, Character>();
		map.put('A', 'T');
		map.put('a', 't');
		map.put('T', 'A');
		map.put('t', 'a');
		map.put('C', 'G');
		map.put('c', 'g');
		map.put('G', 'C');
		map.put('g', 'c');
		
		map.put('N', 'N');
		map.put('n', 'n');
		reverseComplementMap = map;
	}
	public static char reverseComplement(char c){
		return reverseComplementMap.get(c);
	}
}
