package adnascreen;

import java.util.List;

/**
 * Paired read for Reichlab ancient DNA
 * The paired read is a merged result of forward and reverse reads
 * with merging, and barcode and adapter trimming
 * @author Matthew Mah
 *
 */
public class MergedRead extends Read{
	private IndexAndBarcodeKey key;
	
	public static final int maxPassingAlignmentsToConsider = 4;
	public static final int maxQuality = 50;
	public static final char KEY_SEPARATOR = ';';

	public MergedRead(Read r, IndexAndBarcodeKey key){
		super(r);
		this.key = key;
	}
	
	public IndexAndBarcodeKey getKey(){
		return key;
	}
	
	// TODO merged read with barcode and indices in the header should include the capability
	// to convert to String and back again
	/*
	 * To facilitate use with BWA, we include the key string representation with the initial
	 * FASTQ read group. The key representation will be copied to the SAM alignment QNAME field. 
	 * 
	 * For example:
	 * 
	 * TODO
	 * 
	 * Returns a string in FASTQ format with barcode and index key included
	 */
	public String toString(){
		StringBuilder builder = new StringBuilder();
		builder.append(super.toString()); // FASTQ format
		int firstSpaceIndex = builder.indexOf(" ");
		builder.insert(firstSpaceIndex, key.toString());
		builder.insert(firstSpaceIndex, KEY_SEPARATOR);
		return builder.toString();
	}
	
	/**
	 * Find the key of 4-tuple of indices and barcodes for this paired read with index reads.  
	 * Matching is performed using the constraints defined in the argument BarcodeMatchers.
	 * This also checks that the headers of each read match with the exception of the read number
	 * @param r1 forward read
	 * @param r2 reverse read
	 * @param i1 index read for i7
	 * @param i2 index read for i5
	 * @param i5Indices i5 index must be found here to merge
	 * @param i7Indices i7 index must be found here to merge
	 * @param barcodes Inline barcodes must be found here to merge
	 * @return key of 4-tuple of indices and barcodes, or null if all four of these are not found. 
	 */
	public static IndexAndBarcodeKey findExperimentKey(Read r1, Read r2, Read i1, Read i2, 
			BarcodeMatcher i5Indices, BarcodeMatcher i7Indices, BarcodeMatcher barcodes){
		// check for metadata consistency
		if(!r1.getFASTQHeader().equalsExceptRead(r2.getFASTQHeader())
				|| !r1.getFASTQHeader().equalsExceptRead(i1.getFASTQHeader())
				|| !r1.getFASTQHeader().equalsExceptRead(i2.getFASTQHeader()))
			throw new IllegalArgumentException("FASTQ metadata mismatch");

		// Index 1 is i7, Index 2 is i5
		DNASequence i7IndexRaw = i1.getDNASequence();
		DNASequence i5IndexRaw = i2.getDNASequence();
		
		// match indices and barcodes against the known sets to see whether we should process this read pair
		String i5IndexLabel = i5Indices.find(i5IndexRaw);
		String i7IndexLabel = i7Indices.find(i7IndexRaw);

		if(i5IndexLabel != null && i7IndexLabel != null){
			List<Integer> barcodeLengths = barcodes.getBarcodeLengths();
			for(int barcodeLength : barcodeLengths){
				DNASequence p5BarcodeRaw = r1.getDNASequence().subsequence(0, barcodeLength);
				DNASequence p7BarcodeRaw = r2.getDNASequence().subsequence(0, barcodeLength);

				String p5BarcodeLabel = barcodes.find(p5BarcodeRaw);
				String p7BarcodeLabel = barcodes.find(p7BarcodeRaw);

				// process this read pair, which may match an experiment
				if(p5BarcodeLabel != null && p7BarcodeLabel != null){ 
					return new IndexAndBarcodeKey(i5IndexLabel, i7IndexLabel, p5BarcodeLabel, p7BarcodeLabel);
				}
			}
			return new IndexAndBarcodeKey(i5IndexLabel, i7IndexLabel, null, null);
		}
		return null;
	}

	/**
	 * Take a forward and reverse read pair, and merge them if they
	 * have a minimum overlap and resulting length. 
	 * @param r1 forward read
	 * @param r2 reverse read
	 * @param key 4-tuple of indices and barcodes
	 * @param barcodeLength number of base-pairs in a barcode
	 * @param maxPenalty
	 * @param minOverlap Overlap between the forward and reverse reads must be at least this length
	 * @param minMergedLength The resulting merged sequence must be at least this length
	 * @return merged sequence, or null if failure to merge
	 */
	public static MergedRead mergePairedSequences(Read r1, Read r2, IndexAndBarcodeKey key, 
			int r1BarcodeLength, int r2BarcodeLength, int maxPenalty, int minOverlap, int minMergedLength,
			int mismatchPenaltyHigh, int mismatchPenaltyLow, int mismatchBaseQualityThreshold){
		// check for metadata consistency
		if(!r1.getFASTQHeader().equalsExceptRead(r2.getFASTQHeader() ) )
			throw new IllegalArgumentException("FASTQ metadata mismatch");

		// trim leading barcodes and trailing N's
		Read trimmedR1 = r1.subsequence(r1BarcodeLength, r1.length()).trimTrailingUnknownBases();
		Read trimmedR2 = r2.subsequence(r2BarcodeLength, r2.length()).trimTrailingUnknownBases();
		// find best alignment of forward read and reverse-complemented reverse read
		Read reverseComplementR2 = trimmedR2.reverseComplement();
		List<Integer> alignments = Read.findBestAlignment(trimmedR1, reverseComplementR2, maxPenalty, 
				minOverlap, minMergedLength, maxPassingAlignmentsToConsider,
				mismatchPenaltyHigh, mismatchPenaltyLow, mismatchBaseQualityThreshold);
		if(alignments.size() == 1){ // unambiguous merge
			// merge reads
			Read merged = mergeReads(trimmedR1, reverseComplementR2, alignments.get(0));
			MergedRead pairedRead = new MergedRead(merged, key);
			return pairedRead;
		}
		return null;
	}
	
	/**
	 * Merge two reads using the alignment defined by the offset parameter. 
	 * Reads are merged directly; no reverse complement is applied. 
	 * Trimming of adapters is implicit; the merge ends with the 
	 * reverse-read end (beginning of non-reverse complement).  
	 * @param a
	 * @param b
	 * @param offset (a_offset - b_offset)
	 * positive offset:
	 * aaaaaaaaaaaaa
	 *        bbbbbbbbbbbb
	 * negative offset:
	 *        aaaaaaaaaaaa
	 * bbbbbbbbbbbbb
	 * @return
	 */
	static Read mergeReads(Read a, Read b, int offset){
		int resultLength = b.length() + offset;
		StringBuilder dna = new StringBuilder(resultLength);
		int[] quality = new int[resultLength];
		// to loop in three sections
		int i;
		// 1. A only section
		for(i = 0; i < offset; i++){
			dna.append(a.getDNASequence().charAt(i));
			quality[i] = a.getQualitySequence().getQuality(i);
		}
		// 2. merge A, B
		int mergeStart = Math.max(0, offset);
		int mergeEnd = Math.min(a.length(), resultLength);
		for(i = mergeStart; i < mergeEnd; i++){
			char base1 = a.getDNASequence().charAt(i);
			int quality1 = a.getQualitySequence().getQuality(i);
			char base2 = b.getDNASequence().charAt(i - offset);
			int quality2 = b.getQualitySequence().getQuality(i - offset);
			
			BaseWithQuality merged = mergeBases(maxQuality, base1, quality1, base2, quality2);
			dna.append(merged.base);
			quality[i] = merged.quality;
		}
		// 3. B only
		for(i = mergeEnd; i < resultLength; i++){
			dna.append(b.getDNASequence().charAt(i - offset));
			quality[i] = b.getQualitySequence().getQuality(i - offset);
		}
		return new Read(a.getFASTQHeader(), new DNASequence(dna.toString()), new QualitySequence(quality));
	}
	
	
	static class BaseWithQuality{
		char base;
		int quality;
		
		public BaseWithQuality(char base, int quality){
			this.base = base;
			this.quality = quality;
		}
	}
	
	/**
	 * Merge two bases. We assume reads are NOT independent, so the resulting quality will not be better than 
	 * both the input quality values. 
	 * @param maxQuality
	 * @param base1
	 * @param quality1
	 * @param base2
	 * @param quality2
	 * @return
	 */
	static BaseWithQuality mergeBases(int maxQuality, char base1, int quality1, char base2, int quality2){
		char base;
		int quality;
		if(base1 == base2){
			base = Character.toUpperCase(base1);
			quality = Math.min(Math.max(quality1, quality2), maxQuality);
		} else {
			quality = Math.abs(quality2 - quality1);
			base = Character.toUpperCase(quality1 >= quality2 ? base1 : base2);
		}
		return new BaseWithQuality(base, quality);
	}
}
