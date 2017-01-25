package adnascreen;

import static org.junit.Assert.*;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class DNASequenceTests {
	@Rule
	public ExpectedException thrown = ExpectedException.none();
	
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
	public void hashCodeTest(){
		DNASequence s1 = new DNASequence("TGACGCA");
		DNASequence s2 = new DNASequence("TGACGCA");
		assertEquals(s1.hashCode(), s2.hashCode());
	}
	
	@Test
	public void illegalCharacter(){
		thrown.expect(IllegalArgumentException.class);
		DNASequence s1 = new DNASequence("TGACGCX");
		fail("DNASequence contains invalid character that was accepted " + s1.toString()); // should not reach here
	}
}
