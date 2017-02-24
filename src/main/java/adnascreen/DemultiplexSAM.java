package adnascreen;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class DemultiplexSAM {
	public static void main(String [] args) throws IOException {
		Map<IndexAndBarcodeKey, SAMFileWriter> outputFiles = new HashMap<IndexAndBarcodeKey, SAMFileWriter>();
		SAMFileHeader masterHeader = null;
		SAMFileWriterFactory outputFileFactory = new SAMFileWriterFactory();

		// iterate through input files
		for(String filename : args){
			try(
					SamReader reader = SamReaderFactory.makeDefault().open(new File(filename))
			){
				SAMFileHeader header = reader.getFileHeader();
				if(masterHeader == null){
					masterHeader = header;
				} else if(!masterHeader.equals(header)){
					throw new IllegalArgumentException("SAM headers do not match");
				}

				// iterate through alignments
				reader.forEach(record->{
					// parse out key, which is a 4-tuple of indices and barcodes
					String readname = record.getReadName();
					String [] readnameParts = readname.split(String.valueOf(MergedRead.KEY_SEPARATOR));
					IndexAndBarcodeKey key = new IndexAndBarcodeKey(readnameParts[1]);
					// find file corresponding to this key
					SAMFileWriter output = outputFiles.get(key);
					if(output == null){ // open new file, if none exists for this key
						String outputFilename = key.toString() + ".sam";
						File outputFile = new File(outputFilename);
						output = outputFileFactory.makeSAMWriter(header, true, outputFile);
						outputFiles.put(key, output); // 
					}
					// write out to alignment to file
					output.addAlignment(record);
				});
			}
		}
		// cleanup, close all files
		outputFiles.forEach((key, output)->{
			output.close();
		});
	}
}
