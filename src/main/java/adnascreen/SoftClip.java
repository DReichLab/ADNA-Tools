package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

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
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
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
public class SoftClip {
	public static void main(String [] args) throws ParseException, IOException{
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("n", "numBases", true, "Number of bases to clip from both start and end");
		options.addRequiredOption("i", "input BAM", true, "Input BAM filename");
		options.addRequiredOption("o", "output BAM", true, "Output BAM filename");
		options.addOption("b", "BAM", false, "Use bam files for output");
		CommandLine commandLine	= parser.parse(options, args);
		
		int numberOfBasesToClip = Integer.valueOf(commandLine.getOptionValue('n'));
		String inputFilename = commandLine.getOptionValue('i');
		String outputFilename = commandLine.getOptionValue('o');
		boolean useBAM = commandLine.hasOption('b') || Driver.isBAMFilename(outputFilename);
		
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
					softClipBothEndsOfRead(record, numberOfBasesToClip);
					output.addAlignment(record);
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
		int startingAlignmentPosition = record.getAlignmentStart();
		Cigar startingCigar = record.getCigar();
		char [] cigarUnrolled = CigarUtil.cigarArrayFromElements(startingCigar.getCigarElements() );
		numberOfBasesToClip = Math.min(numberOfBasesToClip, cigarUnrolled.length);
		
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
				cigarUnrolled[n] = (char) CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP);
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
				cigarUnrolled[cigarUnrolled.length - n - 1] = (char) CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP);
				clippedBases++;
			}
			else if (!endSide.equals(CigarOperator.SOFT_CLIP)){
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
