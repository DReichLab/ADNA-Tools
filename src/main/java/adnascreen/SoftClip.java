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
		boolean useBAM = commandLine.hasOption('b');
		
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
	 * Soft clip the requested number bases both at the beginning and end of a read.  
	 * If a base is already soft clipped, it still counts towards the requested bases.
	 * Alignment position is adjusted accordingly, taking into account whether the 
	 * cigar operators in the clipped positions consume reference bases. 
	 * @param record
	 * @param numberOfBasesToClip
	 */
	public static void softClipBothEndsOfRead(SAMRecord record, int numberOfBasesToClip){
		int startingAlignmentPosition = record.getAlignmentStart();
		Cigar startingCigar = record.getCigar();
		char [] cigarUnrolled = CigarUtil.cigarArrayFromElements(startingCigar.getCigarElements() );
		numberOfBasesToClip = Math.min(numberOfBasesToClip, cigarUnrolled.length);
		
		int alignmentStartOffset = 0; // number of bases to move alignment start
		for (int n = 0; n < numberOfBasesToClip; n++){
			CigarOperator startSide = CigarOperator.characterToEnum(cigarUnrolled[n]);
			if(startSide.consumesReferenceBases()){
				alignmentStartOffset++;
			}
			cigarUnrolled[n] = (char) CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP);
			cigarUnrolled[cigarUnrolled.length - n - 1] = (char) CigarOperator.enumToCharacter(CigarOperator.SOFT_CLIP);
		}
		
		record.setAlignmentStart(startingAlignmentPosition + alignmentStartOffset);
		String softClippedCigarString = CigarUtil.cigarStringFromArray(cigarUnrolled);
		Cigar softClippedCigar = TextCigarCodec.decode(softClippedCigarString);
		record.setCigar(softClippedCigar);
		// remove MD tag, which is not well documented and needs to match cigar
		record.setAttribute("MD", null);
	}
}
