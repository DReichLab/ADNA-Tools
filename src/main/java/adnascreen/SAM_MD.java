package adnascreen;

import java.util.LinkedList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * For manipulating the optional SAM field MD 
 *
 */
public class SAM_MD {
	// valid MD fields match this regular expression
	public static String VALIDATION = "[0-9]+(([A-Z]|\\^[A-Z]+)[0-9]+)*";
	public static String DIGITS = "[0-9]+";
	public static String MISMATCH = "([A-Z]|\\^[A-Z]+)";
	
	public static Pattern VALIDATION_REGEX = Pattern.compile(VALIDATION);
	public static Pattern DIGITS_REGEX = Pattern.compile(DIGITS);
	public static Pattern MISMATCH_REGEX = Pattern.compile(MISMATCH);
	
	// A, T, C, G reference was this base and mismatched read
	// read inserts will not appear in MD string
	// a, t, c, g reference was this base and was deleted in read
	private LinkedList<Character> positions;
	
	public SAM_MD(String md_field){
		Matcher validator = VALIDATION_REGEX.matcher(md_field);
		if(!validator.matches())
			throw new IllegalArgumentException();
		
		positions = new LinkedList<Character>();
		
		String remaining = md_field.trim();
		
		Matcher digits = null;
		int numMatches = 0;
		
		digits = DIGITS_REGEX.matcher(remaining);
		digits.find();
		numMatches = Integer.valueOf(digits.group());
		remaining = remaining.substring(digits.end());
		for(int n = 0; n < numMatches; n++){
			positions.add('m');
		}
		
		while(remaining.length() > 0){
			Matcher mismatchSectionMatcher = MISMATCH_REGEX.matcher(remaining);
			mismatchSectionMatcher.find();
			String mismatchSection = mismatchSectionMatcher.group();
			remaining = remaining.substring(mismatchSectionMatcher.end());
			if(mismatchSection.charAt(0) == '^'){
				for(int n = 1; n < mismatchSection.length(); n++){ // skip the ^
					positions.add(Character.toLowerCase(mismatchSection.charAt(n)));
				}
			}
			else{
				positions.add(mismatchSection.charAt(0));
				if(mismatchSection.length() > 1)
					throw new IllegalArgumentException("More than mismatch in a row");
			}
			
			digits = DIGITS_REGEX.matcher(remaining);
			digits.find();
			numMatches = Integer.valueOf(digits.group());
			remaining = remaining.substring(digits.end());
			for(int n = 0; n < numMatches; n++){
				positions.add('m');
			}
		}
	}
	
	public void clip(int numLeftElements, int numRightElements){
		// remove right elements from end
		for(int n = 0; n < numRightElements; n++){
			positions.remove(positions.size() - 1);
		}
		// remove left elements from start
		for(int n = 0; n < numLeftElements; n++){
			positions.remove(0);
		}
	}
	
	public int editDistance(){
		int editCount = 0;
		for (char current : positions){
			if(current != 'm'){
				editCount++;
			}
		}
		return editCount;
	}
	
	public String toString(){
		StringBuilder builder = new StringBuilder();
		int matchCount = 0;
		boolean deleteMarker = false;
		for(char current : positions){
			switch(current){
			case 'm':
				matchCount++;
				deleteMarker = false;
				break;
			case 'A':
			case 'C':
			case 'T':
			case 'G':
			case 'N':
				builder.append(matchCount);
				matchCount = 0;
				builder.append(current);
				deleteMarker = false;
				break;
			case 'a':
			case 'c':
			case 't':
			case 'g':
			case 'n':
				if(!deleteMarker){
					builder.append(matchCount);
					builder.append('^');
				}
				matchCount = 0;
				deleteMarker = true;
				builder.append(Character.toUpperCase(current));
				break;
			default:
				throw new IllegalArgumentException("Unexpected MD character " + current);	
			}
		}
		builder.append(matchCount);
		return builder.toString();
	}
}
