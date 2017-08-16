package adnascreen;

import java.util.LinkedList;
import java.util.List;

import htsjdk.samtools.fastq.FastqRecord;

/**
 * 
 * @author Matthew Mah
 *
 */
public class Read {
	private FASTQHeader header; // should probably not be a FASTQHeader, but something more generic
	private DNASequence dna;
	private QualitySequence quality;
	
	public Read(String header, String dna, String quality){
		if(header != null && header.length() > 0 && header.charAt(0) == '@')
			this.header = new FASTQHeader(header);
		else
			this.header = null;
		this.dna = new DNASequence(dna);
		this.quality = new QualitySequence(quality);
		if(dna.length() != quality.length())
			throw new IllegalArgumentException();
	}
	
	public Read(FASTQHeader header, DNASequence dna, QualitySequence quality){
		if(dna == null || quality == null || dna.length() != quality.length())
			throw new IllegalArgumentException();
		this.header = header;
		this.dna = dna;
		this.quality = quality;
		if(dna.length() != quality.length())
			throw new IllegalArgumentException();
	}
	
	public Read(Read other){
		this.header = other.header;
		this.dna = other.dna;
		this.quality = other.quality;
	}
	
	public Read(FastqRecord record){
		this.header = new FASTQHeader(record.getReadName());
		this.dna = new DNASequence(record.getReadString());
		this.quality = new QualitySequence(record.getBaseQualityString());
	}
	
	public int length(){
		return dna.length();
	}
	
	public FASTQHeader getFASTQHeader(){
		return header;
	}
	
	public DNASequence getDNASequence(){
		return dna;
	}
	
	public QualitySequence getQualitySequence(){
		return quality;
	}
	
	// TODO there should be a way to mark a Read as being in reverse order
	public Read reverseComplement(){
		DNASequence reverseComplementDNA = dna.reverseComplement();
		QualitySequence reversedQuality = quality.reverse();
		return new Read(header, reverseComplementDNA, reversedQuality);
	}
	
	public Read subsequence(int beginIndex, int endIndex){
		DNASequence dnaSubsequence = dna.subsequence(beginIndex, endIndex);
		QualitySequence qualitySubsequence = quality.subsequence(beginIndex, endIndex);
		return new Read(header, dnaSubsequence, qualitySubsequence);
	}
	
	public String toString(){
		StringBuilder builder = new StringBuilder();
		builder.append(header);
		builder.append('\n');
		builder.append(dna);
		builder.append('\n');
		builder.append('+');
		builder.append('\n');
		builder.append(quality);
		return builder.toString();
	}
	
	public static Read merge(Read a, Read b, int min_overlap){
		// iterate overlap
		return null;
	}
	
	/**
	 * 
	 * @return
	 */
	public Read trimTrailingUnknownBases(){
		int i = dna.length();
		while((i > 0) && ('N' == Character.toUpperCase(dna.charAt(i-1))) ){
			i--;
		}
		DNASequence trimmedDNA = dna.subsequence(0, i);
		QualitySequence trimmedQuality = quality.subsequence(0, i);
		Read trimmedRead = new Read(header, trimmedDNA, trimmedQuality);
		return trimmedRead;
	}
	
	/**
	 * Evaluate the alignment of two read DNA sequences with offsets
	 * Any reverse complement operation should be called before using this method. 
	 * @param a
	 * @param b
	 * @param aOffset
	 * @param bOffset
	 * @param minOverlapLength
	 * @param maxPenalty defines the max penalty density with minLength
	 * @return true if at or below max penalty density, false otherwise
	 */
	public static boolean alignmentAssessment(Read a, Read b, int aOffset, int bOffset, int minOverlapLength, int maxPenalty,
			int mismatchPenaltyHigh, int mismatchPenaltyLow, int mismatchBaseQualityThreshold){
		int penalty = 0;
		for(int i = 0; (i + aOffset) < a.length() && (i + bOffset) < b.length(); i++){
			int a_i = i + aOffset;
			int b_i = i + bOffset;
			
			char aBase = a.dna.charAt(a_i);
			char bBase = b.dna.charAt(b_i);
			if(aBase != bBase){
				boolean highConfidence = a.quality.getQuality(a_i) >= mismatchBaseQualityThreshold 
						&& b.quality.getQuality(b_i) >= mismatchBaseQualityThreshold;
				penalty += highConfidence ? mismatchPenaltyHigh : mismatchPenaltyLow;
				// cutoff
				if(i <= minOverlapLength && penalty > maxPenalty) return false;
				if(i > minOverlapLength && penalty * minOverlapLength > i * maxPenalty) return false;
			}
		}
		return true;
	}
	
	public static List<Integer> findBestAlignment(Read a, Read b, int maxPenalty, int minOverlapLength, 
			int minResultLength, int maxPositions,
			int mismatchPenaltyHigh, int mismatchPenaltyLow, int mismatchBaseQualityThreshold){
		List<Integer> alignments = new LinkedList<Integer>();
		
		// aOffset and bOffset are the positions to start comparisons in their respective sequences
		// offset = aOffset - bOffset
		// positive:
		// aaaaaaaaaa
		//       bbbbbbbbbb
		// negative:
		//     aaaaaaaaaa
		// bbbbbbbbbbb
		// starting position is for minimum required overlap
		// ending position is for minimum result length
		for(int offset = a.length() - minOverlapLength; b.length() + offset > minResultLength; offset--){
			int aOffset, bOffset;
			if(offset >= 0){
				aOffset = offset;
				bOffset = 0;
			} else{
				aOffset = 0;
				bOffset = Math.abs(offset);
			}
			boolean alignmentAtThisOverlap = alignmentAssessment(a, b, aOffset, bOffset, minOverlapLength, maxPenalty,
					mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
			if(alignmentAtThisOverlap){
				alignments.add(offset);
				if(alignments.size() >= maxPositions){
					return alignments;
				}
			}
		}
		return alignments;
	}
}
