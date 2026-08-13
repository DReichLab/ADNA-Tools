package adnatools;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMTag;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class MarkDuplicatesTest {
	@Test
	public void testAverageBaseQuality() {
		SAMRecord record = new SAMRecord(null);
		record.setBaseQualityString("EEE");
		assertEquals(36.0, MarkDuplicates.averageBaseQuality(record), 0.001);
		record.setBaseQualityString("666");
		assertEquals(21.0, MarkDuplicates.averageBaseQuality(record), 0.001);
	}

	@Test
	public void testDuplicateKeySeparatesStrand() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope();
		SAMRecord forward = new SAMRecord(null);
		forward.setReferenceName("1");
		forward.setAlignmentStart(100);
		forward.setReadBases("ACGT".getBytes());
		forward.setReadNegativeStrandFlag(false);
		forward.setAttribute("XD", "TAG1");

		SAMRecord reverse = new SAMRecord(null);
		reverse.setReferenceName("1");
		reverse.setAlignmentStart(100);
		reverse.setReadBases("ACGT".getBytes());
		reverse.setReadNegativeStrandFlag(true);
		reverse.setAttribute("XD", "TAG1");

		assertFalse(MarkDuplicates.duplicateKey(forward, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(reverse, "XD", true, libraryScope)));
	}

	@Test
	public void testDuplicateKeySeparatesDifferentTagValues() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope();
		SAMRecord tag1 = buildRecord("libA", "1", 100, true, "TAG1");
		SAMRecord tag2 = buildRecord("libA", "1", 100, true, "TAG2");

		assertFalse(MarkDuplicates.duplicateKey(tag1, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(tag2, "XD", true, libraryScope)));
	}

	@Test
	public void testDuplicateKeyIgnoresTagValueWhenTagDisabled() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope();
		SAMRecord tag1 = buildRecord("libA", "1", 100, true, "TAG1");
		SAMRecord tag2 = buildRecord("libA", "1", 100, true, "TAG2");

		assertEquals(MarkDuplicates.duplicateKey(tag1, "XD", false, libraryScope),
				MarkDuplicates.duplicateKey(tag2, "XD", false, libraryScope));
	}

	@Test
	public void testMarkDuplicatesIgnoresTagValuesWhenTagDisabled() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates_tag_values.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_tag_disabled", ".sam");

		try {
			MarkDuplicates markDuplicates = new MarkDuplicates(DemultiplexSAM.duplicatesSAMTag, false, false,
					new MarkDuplicates.LibraryDeduplicationScope());
			markDuplicates.markDuplicates(inputFilename, outputPath.toString(), null, 5, false);

			Map<String, Boolean> duplicateFlags = readDuplicateFlags(outputPath.toString());
			assertEquals(2, duplicateFlags.size());
			assertFalse(duplicateFlags.get("read_tag1"));
			assertTrue(duplicateFlags.get("read_tag2"));
		} finally {
			Files.deleteIfExists(outputPath);
		}
	}

	@Test
	public void testNoTagCommandLineOptionDisablesTagInKey() throws Exception {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates_tag_values.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_no_tag_cli", ".sam");

		try {
			MarkDuplicates.main(new String[] { "-i", inputFilename, "-o", outputPath.toString(), "--no-tag" });

			Map<String, Boolean> duplicateFlags = readDuplicateFlags(outputPath.toString());
			assertEquals(2, duplicateFlags.size());
			assertFalse(duplicateFlags.get("read_tag1"));
			assertTrue(duplicateFlags.get("read_tag2"));
		} finally {
			Files.deleteIfExists(outputPath);
		}
	}

	@Test
	public void testDuplicateKeySeparatesLibraries() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope();
		SAMRecord libA = buildRecord("libA", "1", 100, true, "TAG1");
		SAMRecord libB = buildRecord("libB", "1", 100, true, "TAG1");

		assertFalse(MarkDuplicates.duplicateKey(libA, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(libB, "XD", true, libraryScope)));
	}

	@Test
	public void testDuplicateKeySeparatesReferenceContigs() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope();
		SAMRecord contig1 = buildRecord("libA", "1", 100, true, "TAG1");
		SAMRecord contig2 = buildRecord("libA", "2", 100, true, "TAG1");
		SAMRecord mt = buildRecord("libA", "MT", 100, true, "TAG1");

		assertFalse(MarkDuplicates.duplicateKey(contig1, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(contig2, "XD", true, libraryScope)));
		assertFalse(MarkDuplicates.duplicateKey(contig1, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(mt, "XD", true, libraryScope)));
	}

	@Test
	public void testDuplicateKeySeparatesStartPositions() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope();
		SAMRecord start100 = buildRecord("libA", "1", 100, true, "TAG1");
		SAMRecord start200 = buildRecord("libA", "1", 200, true, "TAG1");

		assertFalse(MarkDuplicates.duplicateKey(start100, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(start200, "XD", true, libraryScope)));
	}

	@Test
	public void testMarkDuplicatesAcrossContigsAndPositions() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates_multi_contig.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_multi_contig", ".sam");
		Path summaryPath = Files.createTempFile("mark_duplicates_multi_contig_summary", ".tsv");

		try {
			MarkDuplicates markDuplicates = new MarkDuplicates(DemultiplexSAM.duplicatesSAMTag, true, false,
					new MarkDuplicates.LibraryDeduplicationScope());
			markDuplicates.markDuplicates(inputFilename, outputPath.toString(), null, 5, false);
			markDuplicates.writeSummary(summaryPath.toString());

			Map<String, Boolean> duplicateFlags = readDuplicateFlags(outputPath.toString());
			assertEquals(8, duplicateFlags.size());
			assertFalse(duplicateFlags.get("read_contig1_pos100_a"));
			assertTrue(duplicateFlags.get("read_contig1_pos100_b"));
			assertFalse(duplicateFlags.get("read_contig1_pos200"));
			assertFalse(duplicateFlags.get("read_contig2_pos100_a"));
			assertTrue(duplicateFlags.get("read_contig2_pos100_b"));
			assertFalse(duplicateFlags.get("read_contig2_pos500"));
			assertFalse(duplicateFlags.get("read_mt_pos1000"));
			assertTrue(duplicateFlags.get("read_mt_pos1000_dup"));

			Map<String, String> summary = readSummary(summaryPath.toString());
			assertEquals("5", summary.get("distinct_aligned_reads"));
			assertEquals("3", summary.get("total_duplicate_reads"));
			assertEquals("2", summary.get("histogram_0"));
			assertEquals("3", summary.get("histogram_1"));
		} finally {
			Files.deleteIfExists(outputPath);
			Files.deleteIfExists(summaryPath);
		}
	}

	@Test
	public void testSamePositionOnDifferentContigsAreNotDuplicates() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates_multi_contig.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_same_pos_diff_contig", ".sam");

		try {
			MarkDuplicates markDuplicates = new MarkDuplicates(DemultiplexSAM.duplicatesSAMTag, true, true,
					new MarkDuplicates.LibraryDeduplicationScope());
			markDuplicates.markDuplicates(inputFilename, outputPath.toString(), null, 5, false);

			List<String> readNames = readReadNames(outputPath.toString());
			assertEquals(5, readNames.size());
			assertTrue(readNames.contains("read_contig1_pos100_a"));
			assertTrue(readNames.contains("read_contig1_pos200"));
			assertTrue(readNames.contains("read_contig2_pos100_a"));
			assertTrue(readNames.contains("read_contig2_pos500"));
			assertTrue(readNames.contains("read_mt_pos1000"));
		} finally {
			Files.deleteIfExists(outputPath);
		}
	}

	@Test
	public void testLibraryGroupCombinesLibraries() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope(
				new String[] { "libA,libB" });
		SAMRecord libA = buildRecord("libA", "1", 100, true, "TAG1");
		SAMRecord libB = buildRecord("libB", "1", 100, true, "TAG1");
		SAMRecord libC = buildRecord("libC", "1", 100, true, "TAG1");

		assertEquals(MarkDuplicates.duplicateKey(libA, "XD", true, libraryScope),
				MarkDuplicates.duplicateKey(libB, "XD", true, libraryScope));
		assertFalse(MarkDuplicates.duplicateKey(libA, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(libC, "XD", true, libraryScope)));
		assertFalse(MarkDuplicates.duplicateKey(libB, "XD", true, libraryScope)
				.equals(MarkDuplicates.duplicateKey(libC, "XD", true, libraryScope)));
	}

	@Test
	public void testLibraryGroupAcceptsSpaceDelimitedLibraries() {
		MarkDuplicates.LibraryDeduplicationScope libraryScope = new MarkDuplicates.LibraryDeduplicationScope(
				new String[] { "libA libB" });
		SAMRecord libA = buildRecord("libA", "1", 100, true, "TAG1");
		SAMRecord libB = buildRecord("libB", "1", 100, true, "TAG1");

		assertEquals(MarkDuplicates.duplicateKey(libA, "XD", true, libraryScope),
				MarkDuplicates.duplicateKey(libB, "XD", true, libraryScope));
	}

	@Test
	public void testMarkDuplicatesKeepsHighestQualityRead() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_out", ".sam");
		Path summaryPath = Files.createTempFile("mark_duplicates_summary", ".tsv");

		try {
			MarkDuplicates markDuplicates = new MarkDuplicates(DemultiplexSAM.duplicatesSAMTag, true, false,
					new MarkDuplicates.LibraryDeduplicationScope());
			markDuplicates.markDuplicates(inputFilename, outputPath.toString(), null, 5, false);
			markDuplicates.writeSummary(summaryPath.toString());

			Map<String, Boolean> duplicateFlags = readDuplicateFlags(outputPath.toString());
			assertEquals(7, duplicateFlags.size());
			assertFalse(duplicateFlags.get("read1"));
			assertTrue(duplicateFlags.get("read2"));
			assertTrue(duplicateFlags.get("read3"));
			assertFalse(duplicateFlags.get("read4"));
			assertFalse(duplicateFlags.get("read5"));
			assertTrue(duplicateFlags.get("read6"));
			assertFalse(duplicateFlags.get("read7"));

			Map<String, String> summary = readSummary(summaryPath.toString());
			assertEquals("4", summary.get("distinct_aligned_reads"));
			assertEquals("3", summary.get("total_duplicate_reads"));
			assertEquals("2", summary.get("histogram_0"));
			assertEquals("1", summary.get("histogram_1"));
			assertEquals("1", summary.get("histogram_2"));
		} finally {
			Files.deleteIfExists(outputPath);
			Files.deleteIfExists(summaryPath);
		}
	}

	@Test
	public void testRemoveDuplicates() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_removed", ".sam");

		try {
			MarkDuplicates markDuplicates = new MarkDuplicates(DemultiplexSAM.duplicatesSAMTag, true, true,
					new MarkDuplicates.LibraryDeduplicationScope());
			markDuplicates.markDuplicates(inputFilename, outputPath.toString(), null, 5, false);

			List<String> readNames = readReadNames(outputPath.toString());
			assertEquals(4, readNames.size());
			assertTrue(readNames.contains("read1"));
			assertTrue(readNames.contains("read4"));
			assertTrue(readNames.contains("read5"));
			assertTrue(readNames.contains("read7"));
		} finally {
			Files.deleteIfExists(outputPath);
		}
	}

	@Test
	public void testDifferentLibrariesAreNotDuplicates() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates_library_groups.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_per_library", ".sam");

		try {
			MarkDuplicates markDuplicates = new MarkDuplicates(DemultiplexSAM.duplicatesSAMTag, true, false,
					new MarkDuplicates.LibraryDeduplicationScope());
			markDuplicates.markDuplicates(inputFilename, outputPath.toString(), null, 5, false);

			Map<String, Boolean> duplicateFlags = readDuplicateFlags(outputPath.toString());
			assertEquals(3, duplicateFlags.size());
			assertFalse(duplicateFlags.get("read1"));
			assertFalse(duplicateFlags.get("read2"));
			assertFalse(duplicateFlags.get("read3"));
		} finally {
			Files.deleteIfExists(outputPath);
		}
	}

	@Test
	public void testLibraryGroupMarksCrossLibraryDuplicates() throws IOException {
		ClassLoader classLoader = getClass().getClassLoader();
		String inputFilename = classLoader.getResource("mark_duplicates/mark_duplicates_library_groups.sam").getPath();
		Path outputPath = Files.createTempFile("mark_duplicates_library_group", ".sam");

		try {
			MarkDuplicates markDuplicates = new MarkDuplicates(DemultiplexSAM.duplicatesSAMTag, true, false,
					new MarkDuplicates.LibraryDeduplicationScope(new String[] { "libA,libB" }));
			markDuplicates.markDuplicates(inputFilename, outputPath.toString(), null, 5, false);

			Map<String, Boolean> duplicateFlags = readDuplicateFlags(outputPath.toString());
			assertEquals(3, duplicateFlags.size());
			assertFalse(duplicateFlags.get("read1"));
			assertTrue(duplicateFlags.get("read2"));
			assertFalse(duplicateFlags.get("read3"));
			assertEquals(2, markDuplicates.getDistinctReadGroups());
			assertEquals(1, markDuplicates.getTotalDuplicateReads());
		} finally {
			Files.deleteIfExists(outputPath);
		}
	}

	private SAMRecord buildRecord(String library, String reference, int start, boolean reverse, String xdTag) {
		String readGroupId = library + "_rg";
		SAMFileHeader header = new SAMFileHeader();
		SAMReadGroupRecord readGroup = new SAMReadGroupRecord(readGroupId);
		readGroup.setLibrary(library);
		readGroup.setSample("sample");
		header.addReadGroup(readGroup);

		SAMRecord record = new SAMRecord(header);
		record.setAttribute(SAMTag.RG.name(), readGroupId);
		record.setReferenceName(reference);
		record.setAlignmentStart(start);
		record.setReadBases("ACGTACGT".getBytes());
		record.setReadNegativeStrandFlag(reverse);
		record.setAttribute("XD", xdTag);
		return record;
	}

	private Map<String, Boolean> readDuplicateFlags(String filename) throws IOException {
		Map<String, Boolean> duplicateFlags = new HashMap<String, Boolean>();
		try (SamReader reader = SamReaderFactory.makeDefault().open(new java.io.File(filename));
				SAMRecordIterator iterator = reader.iterator()) {
			while (iterator.hasNext()) {
				SAMRecord record = iterator.next();
				duplicateFlags.put(record.getReadName(), record.getDuplicateReadFlag());
			}
		}
		return duplicateFlags;
	}

	private List<String> readReadNames(String filename) throws IOException {
		List<String> readNames = new ArrayList<String>();
		try (SamReader reader = SamReaderFactory.makeDefault().open(new java.io.File(filename));
				SAMRecordIterator iterator = reader.iterator()) {
			while (iterator.hasNext()) {
				readNames.add(iterator.next().getReadName());
			}
		}
		return readNames;
	}

	private Map<String, String> readSummary(String filename) throws IOException {
		Map<String, String> summary = new HashMap<String, String>();
		try (BufferedReader reader = new BufferedReader(new FileReader(filename))) {
			String line;
			while ((line = reader.readLine()) != null) {
				String[] fields = line.split("\t");
				if (fields.length == 2 && fields[0].equals("distinct_aligned_reads")) {
					summary.put("distinct_aligned_reads", fields[1]);
				} else if (fields.length == 2 && fields[0].equals("total_duplicate_reads")) {
					summary.put("total_duplicate_reads", fields[1]);
				} else if (fields.length == 2 && fields[0].matches("\\d+")) {
					summary.put("histogram_" + fields[0], fields[1]);
				}
			}
		}
		return summary;
	}
}
