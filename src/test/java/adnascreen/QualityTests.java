package adnascreen;

import static org.junit.Assert.*;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.ExpectedException;

public class QualityTests {
	@Rule
	public ExpectedException thrown = ExpectedException.none();
	
	@Test
	public void q33(){
		assertEquals(0, QualitySequence.quality33Score('!'));
		assertEquals(1, QualitySequence.quality33Score('"'));
		assertEquals(2, QualitySequence.quality33Score('#'));
		assertEquals(3, QualitySequence.quality33Score('$'));
		assertEquals(4, QualitySequence.quality33Score('%'));
		assertEquals(5, QualitySequence.quality33Score('&'));
		assertEquals(6, QualitySequence.quality33Score('\''));
		assertEquals(7, QualitySequence.quality33Score('('));
		assertEquals(8, QualitySequence.quality33Score(')'));
		assertEquals(9, QualitySequence.quality33Score('*'));
		assertEquals(10, QualitySequence.quality33Score('+'));
		assertEquals(11, QualitySequence.quality33Score(','));
		assertEquals(12, QualitySequence.quality33Score('-'));
		assertEquals(13, QualitySequence.quality33Score('.'));
		assertEquals(14, QualitySequence.quality33Score('/'));
		assertEquals(15, QualitySequence.quality33Score('0'));
		assertEquals(16, QualitySequence.quality33Score('1'));
		assertEquals(17, QualitySequence.quality33Score('2'));
		assertEquals(18, QualitySequence.quality33Score('3'));
		assertEquals(19, QualitySequence.quality33Score('4'));
		assertEquals(20, QualitySequence.quality33Score('5'));
		assertEquals(21, QualitySequence.quality33Score('6'));
		assertEquals(22, QualitySequence.quality33Score('7'));
		assertEquals(23, QualitySequence.quality33Score('8'));
		assertEquals(24, QualitySequence.quality33Score('9'));
		assertEquals(25, QualitySequence.quality33Score(':'));
		assertEquals(26, QualitySequence.quality33Score(';'));
		assertEquals(27, QualitySequence.quality33Score('<'));
		assertEquals(28, QualitySequence.quality33Score('='));
		assertEquals(29, QualitySequence.quality33Score('>'));
		assertEquals(30, QualitySequence.quality33Score('?'));
		assertEquals(31, QualitySequence.quality33Score('@'));
		assertEquals(32, QualitySequence.quality33Score('A'));
		assertEquals(33, QualitySequence.quality33Score('B'));
		assertEquals(34, QualitySequence.quality33Score('C'));
		assertEquals(35, QualitySequence.quality33Score('D'));
		assertEquals(36, QualitySequence.quality33Score('E'));
		assertEquals(37, QualitySequence.quality33Score('F'));
		assertEquals(38, QualitySequence.quality33Score('G'));
		assertEquals(39, QualitySequence.quality33Score('H'));
		assertEquals(40, QualitySequence.quality33Score('I'));
		
		assertEquals('!', QualitySequence.quality33Char(0));
		assertEquals('"', QualitySequence.quality33Char(1));
	}
	
	@Test
	public void reverse(){
		int [] scores = {1, 2, 3, 4, 5, 6, 7, 8, 9};
		QualitySequence forward = new QualitySequence(scores);
		QualitySequence backward = forward.reverse();
		for(int i = 0; i < forward.length(); i++){
			assertEquals(forward.getQuality(i), backward.getQuality(backward.length() - 1 - i));
			assertEquals(backward.getQuality(i), forward.getQuality(forward.length() - 1 - i));
		}
	}
	
	@Test
	public void qualitySample(){
		String quality = "AE/E<EEAE<E<E";
		QualitySequence q = new QualitySequence(quality);
		String converted = q.toString();
		assertEquals(quality, converted);
	}
	
	@Test
	public void badBounds() {
		thrown.expect(ArrayIndexOutOfBoundsException.class);
		String quality = "AE/E<EEAE<E<E";
		QualitySequence q = new QualitySequence(quality);
		q.subsequence(0, 20);
	}
	
	@Test
	public void isEqual(){
		String quality = "AE/E<EEAE<E<E";
		QualitySequence q1 = new QualitySequence(quality);
		QualitySequence q2 = new QualitySequence(new String(quality));
		assertEquals(q1, q2);
		assertEquals(q1.toString(), q2.toString());
	}
	
	@Test
	public void isNotEqual(){
		String quality1 = "AE/E<EEAE<E<E";
		String quality2 = "BE/E<EEAE<E<E";
		QualitySequence q1 = new QualitySequence(quality1);
		QualitySequence q2 = new QualitySequence(quality2);
		assertNotEquals(q1, q2);
		assertNotEquals(q1.toString(), q2.toString());
	}
}
