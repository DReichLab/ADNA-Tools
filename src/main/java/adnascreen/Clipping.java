package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMValidationError;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.TextCigarCodec;
import htsjdk.samtools.util.CigarUtil;

/**
 * To prevent damaged bases from affecting downstream analysis, 
 * we soft clip bases from the ends of reads
 * 
 * @author mmah
 *
 */
public class Clipping {
	private int defaultNumberOfBasesToClip;
	private HashMap<String, Integer> clippingLengthByLibrary;
	
	/**
	 * Parse command line options for clipping lengths based on library
	 * This is designed for 3 options (UDG minus, half, plus)
	 * Options in the command line need to match those in addSoftClipCommandLineOptions
	 * @param commandLine
	 */
	public Clipping(CommandLine commandLine) {
		defaultNumberOfBasesToClip = Integer.valueOf(commandLine.getOptionValue('n', "0"));
		clippingLengthByLibrary = new HashMap<String, Integer>();
		
		int numberOfBasesToClipSpecial1 = Integer.valueOf(commandLine.getOptionValue('x', "0"));
		String[] specialClippingLibrariesArray1 = commandLine.getOptionValues("libraries1");
		if (specialClippingLibrariesArray1 != null) {
			for(String library : specialClippingLibrariesArray1) {
				clippingLengthByLibrary.put(library, numberOfBasesToClipSpecial1);
			}
		}
		int numberOfBasesToClipSpecial2 = Integer.valueOf(commandLine.getOptionValue('y', "0"));
		String[] specialClippingLibrariesArray2 = commandLine.getOptionValues("libraries2");
		if (specialClippingLibrariesArray2 != null) {
			for(String library : specialClippingLibrariesArray2) {
				clippingLengthByLibrary.put(library, numberOfBasesToClipSpecial2);
			}
		}
	}
	
	/**
	 * 
	 * @param record
	 * @return number of bases to clip on both ends based on library
	 */
	public int getClippingLength(SAMRecord record) {
		int numberOfBasesToClip = defaultNumberOfBasesToClip;
		// clipping length depends on library
		SAMReadGroupRecord readGroup = record.getReadGroup();
		if(readGroup != null) {
			String library = readGroup.getLibrary();
			if(library != null) {
				numberOfBasesToClip = clippingLengthByLibrary.getOrDefault(library, defaultNumberOfBasesToClip);
			}
		}
		return numberOfBasesToClip;
	}
	
	public static void addSoftClipCommandLineOptions(Options options) {
		options.addOption("n", "numBases", true, "Number of bases to clip from both start and end");
		// For mixed UDG bams, specify the number of bases for a second and third case 
		options.addOption("x", "numBasesSpecial1", true, "Number of bases to clip from both start and end for reads in corresponding read groups, as exception to normal clipping");
		options.addOption("s", "libraries1", true, "libraries to clip with numBasesSpecial1");
		options.addOption("y", "numBasesSpecial2", true, "Number of bases to clip from both start and end for reads in corresponding read groups, as exception to normal clipping");
		options.addOption("t", "libraries2", true, "libraries to clip with numBasesSpecial2");
	}
	
	// non-empty reads have at least one (mis)match
	public static boolean isNonEmptyRead(SAMRecord record) {
		Cigar cigar = record.getCigar();
		return cigar.containsOperator(CigarOperator.MATCH_OR_MISMATCH);
	}
	
	public static void main(String [] args) throws ParseException, IOException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input BAM", true, "Input BAM filename");
		options.addRequiredOption("o", "output BAM", true, "Output BAM filename");
		options.addOption("b", "BAM", false, "Use bam files for output");
		options.addOption(null, "hard", false, "Hard clip bases. Use this for contammix, for example, which requires realignment");
		addSoftClipCommandLineOptions(options);
		CommandLine commandLine	= parser.parse(options, args);
		
		String inputFilename = commandLine.getOptionValue('i');
		String outputFilename = commandLine.getOptionValue('o');
		boolean useBAM = commandLine.hasOption('b') || Driver.isBAMFilename(outputFilename);
		boolean hardClip = commandLine.hasOption("hard");
		
		// clipping for command line specified libraries to handle separate UDG treatments 
		Clipping softClipLengths = new Clipping(commandLine);
		
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename)));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
				){
			SAMFileHeader header = reader.getFileHeader();
			
			SAMFileWriter output;
			SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
			BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputFilename));
			if(useBAM){
				output = outputFileFactory.makeBAMWriter(header, false, outputFile);
			} else {
				output = outputFileFactory.makeSAMWriter(header, false, outputFile);
			}

			SAMRecordIterator i = reader.iterator();
			while(i.hasNext()){
				// iterate through alignments
				try{
					SAMRecord record = i.next();
					int numberOfBasesToClip = softClipLengths.getClippingLength(record);				
					if(hardClip)
						hardClipBothEndsOfRead(record, numberOfBasesToClip);
					else
						softClipBothEndsOfRead(record, numberOfBasesToClip);
					
					// very short reads may be entirely clipped, remove these from the data set
					if(isNonEmptyRead(record)) {
						output.addAlignment(record);
						
						// diagnostic prints for validation failures
						List<SAMValidationError> errors = record.isValid();
						if(errors != null) {
							for(SAMValidationError error : errors) {
								System.err.println(error.toString());
							}
							System.err.println(record.toString());
						}
					}
				}
				catch(Exception e){
					System.err.println(e.toString());
				}
			}
			
			output.close();
		}
	}
	
	/**
	 * Edit distance of this cigar from inserts
	 * Other portion of edit distance comes from MD field
	 * @param cigar
	 * @return
	 */
	public static int editDistance(Cigar cigar){
		int edits = 0;
		for(CigarElement c : cigar.getCigarElements()){
			if(c.getOperator().equals(CigarOperator.INSERTION)){
				edits += c.getLength();
			}
		}
		return edits;
	}
	
	public static int editDistance(Cigar cigar, SAM_MD md){
		return editDistance(cigar) + md.editDistance();
	}
	
	/**
	 * Soft clip the requested number bases both at the beginning and end of a read.  
	 * If a base is already soft clipped, it still counts towards the requested bases.
	 * Alignment position is adjusted accordingly, taking into account whether the 
	 * cigar operators in the clipped positions consume reference bases. 
	 * Modify the MD field to match the changed cigar.
	 * Modify the edit distance tag NM to be consistent with the cigar and MD. 
	 * @param record
	 * @param numberOfBasesToClip
	 */
	public static void softClipBothEndsOfRead(SAMRecord record, int numberOfBasesToClip){
		clipBothEndsOfRead(record, numberOfBasesToClip, false);
	}
	
	/**
	 * Soft-clipped bases remain in record, which is what we normally want. 
	 * For reverting back to fastq, however, we want to stop the soft-clipped bases from being reintroduced. 
	 * For this purpose, hard clip bases at the ends of the read. 
	 */
	public static void hardClipBothEndsOfRead(SAMRecord record, int numberOfBasesToClip){
		clipBothEndsOfRead(record, numberOfBasesToClip, true);
	}
		
	public static void clipBothEndsOfRead(SAMRecord record, int numberOfBasesToClip, boolean hardClip){
		int startingAlignmentPosition = record.getAlignmentStart();
		Cigar startingCigar = record.getCigar();
		char [] cigarUnrolled = CigarUtil.cigarArrayFromElements(startingCigar.getCigarElements() );
		numberOfBasesToClip = Math.min(numberOfBasesToClip, cigarUnrolled.length);
		
		CigarOperator replacementCigarOperator = hardClip ? CigarOperator.HARD_CLIP : CigarOperator.SOFT_CLIP;
		
		int alignmentStartOffset = 0; // number of bases to move alignment start, and amount to change MD field
		int endSideChange = 0; // modifies MD only
		
		int n = 0;
		int clippedBases = 0;
		int hardClipped = 0;
		// clip bases from front
		while(clippedBases < numberOfBasesToClip) {
			CigarOperator startSide = CigarOperator.characterToEnum(cigarUnrolled[n]);
			if(startSide.consumesReferenceBases()){
				alignmentStartOffset++;
			}
			if(startSide.consumesReadBases()) {
				cigarUnrolled[n] = (char) CigarOperator.enumToCharacter(replacementCigarOperator);
				if (hardClip)
					hardClipped++;
				clippedBases++;
			}
			else if (!startSide.equals(CigarOperator.SOFT_CLIP)){
				cigarUnrolled[n] = (char) CigarOperator.enumToCharacter(CigarOperator.HARD_CLIP);
				hardClipped++;
			}
			n++;
		}
		// hard clipped bases must be before soft clipped bases
		int totalClippedBases = n;
		for(n = 0; n < totalClippedBases; n++) {
			if(n < hardClipped)
				cigarUnrolled[n] = (char) CigarOperator.enumToCharacter(CigarOperator.HARD_CLIP);
			else
				cigarUnrolled[n] = (char) CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP);
		}
		
		// clip bases from rear
		n = 0;
		clippedBases = 0;
		hardClipped = 0;
		while(clippedBases < numberOfBasesToClip) {
			CigarOperator endSide = CigarOperator.characterToEnum(cigarUnrolled[cigarUnrolled.length - n - 1]);
			if(endSide.consumesReferenceBases()){
				endSideChange++;
			}
			if(endSide.consumesReadBases()) {
				cigarUnrolled[cigarUnrolled.length - n - 1] = (char) CigarOperator.enumToCharacter(replacementCigarOperator);
				if (hardClip)
					hardClipped++;
				clippedBases++;
			}
			else {
				cigarUnrolled[cigarUnrolled.length - n - 1] = (char) CigarOperator.enumToCharacter(CigarOperator.HARD_CLIP);
				hardClipped++;
			}
			n++;
		}
		// hard clipped bases must be before soft clipped bases
		totalClippedBases = n;
		for(n = 0; n < totalClippedBases; n++) {
			if(n < hardClipped)
				cigarUnrolled[cigarUnrolled.length - n - 1] = (char) CigarOperator.enumToCharacter(CigarOperator.HARD_CLIP);
			else
				cigarUnrolled[cigarUnrolled.length - n - 1] = (char) CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP);
		}
		
		// remove hard clipped query bases from the base and quality strings
		if(hardClip) {
			String qualities = record.getBaseQualityString();
			String shortenedQualities = qualities.substring(numberOfBasesToClip, qualities.length() - numberOfBasesToClip);
			record.setBaseQualityString(shortenedQualities);
			String bases = record.getReadString();
			String shortenedBases = bases.substring(numberOfBasesToClip, bases.length() - numberOfBasesToClip);
			record.setReadString(shortenedBases);
		}
		
		record.setAlignmentStart(startingAlignmentPosition + alignmentStartOffset);
		String softClippedCigarString = CigarUtil.cigarStringFromArray(cigarUnrolled);
		Cigar softClippedCigar = TextCigarCodec.decode(softClippedCigarString);
		record.setCigar(softClippedCigar);
		
		// modify the MD field to match cigar
		String mdString = (String) record.getAttribute("MD");
		SAM_MD md = new SAM_MD(mdString);
		md.clip(alignmentStartOffset, endSideChange);
		record.setAttribute("MD", md.toString());
		
		// modify NM edit distance field
		int edit = editDistance(softClippedCigar, md);
		record.setAttribute("NM", edit);
	}
}
