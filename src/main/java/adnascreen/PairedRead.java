package adnascreen;

import java.util.List;

/**
 * Paired read for Reichlab ancient DNA
 * The paired read is a merged result of forward and reverse reads
 * with merging, and barcode and adapter trimming
 * @author Matthew Mah
 *
 */
public class PairedRead extends Read{
	private IndexAndBarcodeKey key;
	
	public static final int maxPassingAlignmentsToConsider = 4;
	public static final int maxQuality = 50;

	public PairedRead(Read r, IndexAndBarcodeKey key){
		super(r);
		this.key = key;
	}
	
	public IndexAndBarcodeKey getKey(){
		return key;
	}
	
	/**
	 * Returns a string in FASTQ format with barcode and index key included
	 */
	public String toString(){
		StringBuilder builder = new StringBuilder();
		builder.append(super.toString());
		int firstNewlineIndex = builder.indexOf("\n");
		builder.insert(firstNewlineIndex, key.toString());
		builder.insert(firstNewlineIndex, " ");
		return builder.toString();
	}
	
	/**
	 * Take a forward and reverse read pair, and 
	 * @param r1
	 * @param r2
	 * @param i1
	 * @param i2
	 * @param i5Indices
	 * @param i7Indices
	 * @param barcodes
	 * @param minOverlap
	 * @param minMergedLength The resulting merged sequence must be 
	 * @return merged sequence, or null if failure to merge
	 */
	public static PairedRead mergePairedSequences(Read r1, Read r2, Read i1, Read i2, 
			IndexMatcher i5Indices, IndexMatcher i7Indices, BarcodeMatcher barcodes, 
			int maxPenalty, int minOverlap, int minMergedLength){
		// check for metadata consistency
		if(!r1.getFASTQHeader().equalsExceptRead(r2.getFASTQHeader())
				|| !r1.getFASTQHeader().equalsExceptRead(i1.getFASTQHeader())
				|| !r1.getFASTQHeader().equalsExceptRead(i2.getFASTQHeader()))
			throw new IllegalArgumentException("FASTQ metadata mismatch");
		
		// Index 1 is i7, Index 2 is i5
		DNASequence i7IndexRaw = i1.getDNASequence();
		DNASequence i5IndexRaw = i2.getDNASequence();
		
		DNASequence p5BarcodeRaw = r1.getDNASequence().subsequence(0, 7);
		DNASequence p7BarcodeRaw = r2.getDNASequence().subsequence(0, 7);
		
		// match indices and barcodes against the known sets to see whether we should process this read pair
		DNASequence i5Index = i5Indices.find(i5IndexRaw);
		DNASequence i7Index = i7Indices.find(i7IndexRaw);
		String p5BarcodeSet = barcodes.find(p5BarcodeRaw);
		String p7BarcodeSet = barcodes.find(p7BarcodeRaw);
		
		if(i5Index != null 
				&& i7Index != null
				&& p5BarcodeSet != null
				&& p7BarcodeSet != null){ // process this read pair, which may match an experiment
			// trim leading barcodes and trailing N's
			Read trimmedR1 = r1.subsequence(barcodes.getBarcodeLength(), r1.length()).trimTrailingUnknownBases();
			Read trimmedR2 = r2.subsequence(barcodes.getBarcodeLength(), r2.length()).trimTrailingUnknownBases();
			// find best alignment of forward read and reverse-complemented reverse read
			Read reverseComplementR2 = trimmedR2.reverseComplement();
			List<Integer> alignments = Read.findBestAlignment(trimmedR1, reverseComplementR2, maxPenalty, 
					minOverlap, minMergedLength, maxPassingAlignmentsToConsider);
			if(alignments.size() == 1){ // unambiguous merge
				// merge reads
				Read merged = mergeReads(trimmedR1, reverseComplementR2, alignments.get(0));
				IndexAndBarcodeKey key = new IndexAndBarcodeKey(i5Index, i7Index, p5BarcodeSet, p7BarcodeSet);
				PairedRead pairedRead = new PairedRead(merged, key);
				return pairedRead;
			}
		}
		return null;
	}
	
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
