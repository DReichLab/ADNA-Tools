package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Set;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.util.AbstractRecordAndOffset;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.tribble.bed.BEDCodec;
import htsjdk.tribble.bed.BEDFeature;

// Filter for 1240k target set with minimum base and mapping qualities
// Allow soft clipping in this pass
public class FilterSAM {
	OverlapDetector<BEDFeature> features;
	
	public static void main(String [] args) throws ParseException, IOException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("i", "input_BAM", true, "Input BAM filename");
		options.addRequiredOption("o", "output_BAM", true, "Output BAM filename");
		options.addRequiredOption("m", "minimum_mapping_quality", true, "minimum mapping quality");
		options.addRequiredOption("q", "minimum_base_quality", true, "minimum base quality");
		options.addRequiredOption("p", "positions", true, "positions file in BED format");
		options.addOption("b", "BAM", false, "Use bam files for output");
		//options.addOption("c", "soft_clip", true, "bases to soft clip from clip for deamination damage");
		SoftClip.addSoftClipCommandLineOptions(options);
		
		CommandLine commandLine	= parser.parse(options, args);
		
		int minimumMappingQuality = Integer.valueOf(commandLine.getOptionValue("minimum_mapping_quality"));
		int minimumBaseQuality = Integer.valueOf(commandLine.getOptionValue("minimum_base_quality"));
		String inputFilename = commandLine.getOptionValue('i');
		String outputFilename = commandLine.getOptionValue('o');
		String bedFilename = commandLine.getOptionValue('p');
		boolean useBAM = commandLine.hasOption('b') || Driver.isBAMFilename(outputFilename);
		SoftClip softClipLengths = new SoftClip(commandLine);
		
		SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(inputFilename)));
		try(
				SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
				){
			FilterSAM filter = new FilterSAM(bedFilename);
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
				SAMRecord record = null;
				try{
					record = i.next();
					if(!record.getReadUnmappedFlag()) {
						int softClipBases = softClipLengths.getClippingLength(record);
						if(softClipBases > 0)
							SoftClip.softClipBothEndsOfRead(record, softClipBases);
						if (SoftClip.isNonEmptyRead(record) && filter.filter(record, minimumMappingQuality, minimumBaseQuality)) {
							output.addAlignment(record);
						}
					}				
				}
				catch(Exception e){
					System.err.println(e.toString());
					if(record != null) {
						System.err.println(record.toString());
					}
				}
			}
			
			output.close();
		}
		finally{}
	}
	
	public FilterSAM(String filename) throws IOException {
		this.readBEDFile(filename);
	}
	
	public void readBEDFile(String filename) throws IOException{
		features = new OverlapDetector<BEDFeature>(0, 0);
		BEDCodec bedCodec = new BEDCodec();
		try(InputStream referenceStream = new FileInputStream(filename);
				BufferedReader reader = new BufferedReader(new InputStreamReader(referenceStream));
				){
			String line;
			while((line = reader.readLine()) != null){
				BEDFeature feature = bedCodec.decode(line);
				features.addLhs(feature, feature);
			}
		}
	}
	
	/**
	 * 
	 * @param record
	 * @param minimumMappingQuality
	 * @param minimumBaseQuality
	 * @return true if record meets quality requirements at a feature
	 */
	public boolean filter(SAMRecord record, int minimumMappingQuality, int minimumBaseQuality) {
		if(record.getMappingQuality() >= minimumMappingQuality) {
			Set<BEDFeature> overlappingFeatures = features.getOverlaps(record);
			for(BEDFeature feature : overlappingFeatures) {
				boolean sufficientQuality = true;
				for(int referencePosition = feature.getStart(); referencePosition <= feature.getEnd(); referencePosition++) {
					int readPosition = record.getReadPositionAtReferencePosition(referencePosition);
					if(readPosition <= 0) { // 1-based, 0 indicates deletion
						sufficientQuality = false;
					} else {
						AbstractRecordAndOffset baseToCheck = new AbstractRecordAndOffset(record, readPosition-1); // offset is 0-based
						if(baseToCheck.getBaseQuality() < minimumBaseQuality)
							sufficientQuality = false;
					}
				}
				if(sufficientQuality)
					return true;
			}
		}
		return false;
	}
}
