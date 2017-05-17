package adnascreen;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMFormatException;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class DemultiplexSAM {
	public static final String ALIGNED = "aligned";
	public static final String DEMULTIPLEXED = "demultiplexed";
	
	public static void main(String [] args) throws IOException, ParseException {
		CommandLineParser parser = new DefaultParser();
		Options options = new Options();
		options.addRequiredOption("s", "statisticsFilename", true, "Statistics file sorted in order of output");
		options.addOption("n", "numSamples", true, "Number of top samples to output");
		options.addOption("b", "BAM", false, "Use bam files for output");
		CommandLine commandLine	= parser.parse(options, args);
		
		int numTopSamples = Integer.valueOf(commandLine.getOptionValue('n', "1000"));
		boolean useBAM = commandLine.hasOption('b');
		String fileExtension = useBAM ? ".bam" : ".sam";
		
		String duplicatesSAMTag = "XD";
		Map<IndexAndBarcodeKey, SAMFileWriter> outputFiles = new HashMap<IndexAndBarcodeKey, SAMFileWriter>(numTopSamples);
		
		SAMSequenceDictionary alignmentReference = null;
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();
		
		// read statistics file with top keys
		// open a SAM/BAM file for each key for demultiplexing its data
		String statisticsFilename = commandLine.getOptionValue('s');
		File f = new File(statisticsFilename);
		try(BufferedReader reader = new BufferedReader(new FileReader(f))){
			Integer.valueOf(reader.readLine());
			for(int n = 0; n < numTopSamples; n++){
				String entryLine = reader.readLine();
				String [] fields = entryLine.split("\t");
				String keyString = fields[0];
				IndexAndBarcodeKey key = new IndexAndBarcodeKey(keyString);
				outputFiles.put(key, null); // mark this key for output later
				// we delay opening SAM/BAM file writer until the SAM/BAM header is available
				// this is after we have opened the first SAM/BAM input file
			}
		}
		SampleSetsCounter statistics = new SampleSetsCounter(f);
		// iterate through input files
		
		List<String> samFilenamesToProcess = commandLine.getArgList();
		for(String filename : samFilenamesToProcess){
			SamInputResource bufferedSAMFile = SamInputResource.of(new BufferedInputStream(new FileInputStream(filename)));
			try(
					SamReader reader = SamReaderFactory.makeDefault().open(bufferedSAMFile);
			){
				SAMFileHeader header = reader.getFileHeader();
				SAMSequenceDictionary currentAlignmentReference = header.getSequenceDictionary();
				if(alignmentReference == null){
					alignmentReference = currentAlignmentReference;
				} else if(!alignmentReference.equals(currentAlignmentReference)){
					throw new IllegalArgumentException("SAM references do not match");
				}

				SAMRecordIterator i = reader.iterator();
				while(i.hasNext()){
					// iterate through alignments
					try{
						SAMRecord record = i.next();
						// parse out key, which is a 4-tuple of indices and barcodes
						String readName = record.getReadName();
						String [] readNameParts = readName.split(String.valueOf(MergedRead.KEY_SEPARATOR));
						IndexAndBarcodeKey key = new IndexAndBarcodeKey(readNameParts[1]);
						IndexAndBarcodeKey keyFlattened = key.flatten();
						
						// for deduplication, write unflattened key identifying barcodes within set and read length to tag
						int length = record.getReadLength();
						record.setAttribute(duplicatesSAMTag, key.toString() + "_" + length);
						// remove the key from the read name
						record.setReadName(readNameParts[0]);

						// record statistics
						// count of demultiplexed reads is for checking consistency
						statistics.increment(keyFlattened, DEMULTIPLEXED);
						if(!record.getReadUnmappedFlag()){ // read is mapped
							statistics.increment(keyFlattened, ALIGNED);
						}
						
						// write only to open files for top keys
						if(outputFiles.containsKey(keyFlattened)){
							// find file corresponding to this key
							SAMFileWriter output = outputFiles.get(keyFlattened);
							if(output == null){ // open new file, if none exists for this key
								String outputFilename = keyFlattened.toString() + fileExtension;
								BufferedOutputStream outputFile = new BufferedOutputStream(new FileOutputStream(outputFilename));
								if(useBAM){
									output = outputFileFactory.makeBAMWriter(header, false, outputFile);
								} else {
									output = outputFileFactory.makeSAMWriter(header, false, outputFile);
								}
								outputFiles.put(keyFlattened, output); // 
							}
							// write alignment to file
							output.addAlignment(record);
						}
					} catch (SAMFormatException e){
						System.err.print(filename + "\t");
						System.err.println(e);
						// ignore this record and continue to the next
					} catch (Exception e){
						System.err.print(filename + "\t");
						System.err.println(e);
					}
				}
			}
		}
		System.out.println(statistics.toStringSorted(IndexAndBarcodeScreener.RAW));
		// cleanup, close all files
		for(SAMFileWriter writer : outputFiles.values()){
			if(writer != null){
				writer.close();
			}
		}
	}
}
