package adnascreen;

import static org.junit.Assert.*;

import org.junit.Test;

public class DNASequenceTests {
	@Test
	public void reverseComplement_ACTGN(){
		String input = "ACTGN";
		String expected = "NCAGT";
		
		DNASequence A = new DNASequence(input);
		DNASequence reverseComplement = A.reverseComplement();
		DNASequence rereverseComplement = reverseComplement.reverseComplement();
		
		assertEquals(reverseComplement.toString(), expected);
		assertEquals(A, rereverseComplement);
		assertNotEquals(A, input); // DNASequence v. String
	}
	
	@Test
	public void exampleReverseComplement(){
		DNASequence A = new DNASequence("ACTGCGTCTAGCATTACTTATATGATATGTCTCCATACCAATTACAATCTCCAAGTGAACGAGATCGGAAGAGCAC");
		DNASequence r = A.reverseComplement();
		DNASequence expected = new DNASequence("GTGCTCTTCCGATCTCGTTCACTTGGAGATTGTAATTGGTATGGAGACATATCATATAAGTAATGCTAGACGCAGT");
		assertEquals(expected, r);
	}
	
	@Test
	public void hamming1(){
		DNASequence s1 = new DNASequence("TGACGCA");
		DNASequence s2 = new DNASequence("AGACGCA");
		assertEquals(s1.hammingDistance(s2), 1);
	}
	
	@Test
	public void hamming7(){
		DNASequence s1 = new DNASequence("TGACGCA");
		DNASequence s2 = new DNASequence("ATCGTGC");
		assertEquals(s1.hammingDistance(s2), 7);
	}
	
	@Test
	public void hammingLengthMismatch() {
		DNASequence s1 = new DNASequence("TGACGCA");
		DNASequence s2 = new DNASequence("ATCGTG");
		assertThrows(IllegalArgumentException.class, () -> s1.hammingDistance(s2));
	}
	
	@Test
	public void hashCodeTest(){
		DNASequence s1 = new DNASequence("TGACGCA");
		DNASequence s2 = new DNASequence("TGACGCA");
		assertEquals(s1.hashCode(), s2.hashCode());
	}
	
	@Test
	public void illegalCharacter(){
		assertThrows(IllegalArgumentException.class, () -> new DNASequence("TGACGCX"));
	}
}
