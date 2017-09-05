package adnascreen;

/**
 * This class allows demultiplexing of experiments. 
 * The IndexAndBarcodeKey comprises a p5 index (i5), a p7 index (i7), a p5 barcode (p5), and a p7 barcode (p7)
 * The p5 and p7 barcodes may optionally be null. 
 */
public class IndexAndBarcodeKey {
	// labels
	private String i5, i7;
	private String p5, p7;

	// _ is chosen because it is a natural character in filenames
	// and it is not contained in any of the index or barcode names  
	public static final char FIELD_SEPARATOR = '_'; // this needs to be distinct from MergedRead.KEY_SEPARATOR
	private static String INDEX_SEPARATOR;
	private static String INDEX_SEPARATOR_REGEX;
	private static final String FIELD_SEPARATOR_STRING = String.valueOf(FIELD_SEPARATOR);
	static{
		INDEX_SEPARATOR = String.valueOf(BarcodeMatcher.INDEX_DELIMITER);
		INDEX_SEPARATOR_REGEX = "\\" + INDEX_SEPARATOR;
	}
	
	public IndexAndBarcodeKey(String i5, String i7, String p5, String p7){
		if(i5 != null && i5.contains(FIELD_SEPARATOR_STRING)
				|| i7 != null && i7.contains(FIELD_SEPARATOR_STRING)
				|| p5 != null && p5.contains(FIELD_SEPARATOR_STRING)
				|| p7 != null && p7.contains(FIELD_SEPARATOR_STRING)){
			throw new IllegalArgumentException("Index or barcode contains illegal separator character " + FIELD_SEPARATOR);
		}
		this.i5 = i5;
		this.i7 = i7;
		this.p5 = p5;
		this.p7 = p7;
	}
	
	public IndexAndBarcodeKey(String keyString){
		String regexPattern = String.valueOf(FIELD_SEPARATOR);
		String [] fields = keyString.split(regexPattern);
		i5 = fields[0];
		i7 = fields[1];
		// barcodes may be null
		// split omits trailing empty strings
		p5 = (fields.length >= 3 && fields[2].length() > 0) ? fields[2] : null;
		p7 = (fields.length >= 4 && fields[3].length() > 0) ? fields[3] : null;
	}
	
	@Override
	public String toString(){
		StringBuilder b = new StringBuilder();
		b.append(i5);
		b.append(FIELD_SEPARATOR);
		b.append(i7);
		b.append(FIELD_SEPARATOR);
		b.append(p5 != null ? p5 : "");
		b.append(FIELD_SEPARATOR);
		b.append(p7 != null ? p7 : "");
		return b.toString();
	}
	
	/**
	 * For each index and barcode, remove the portion of the label that specifies the position
	 * within the reference barcode set.   
	 * This information is necessary for deduplication, but should be removed for demultiplexing. 
	 * @return
	 */
	public IndexAndBarcodeKey flatten(){
		return new IndexAndBarcodeKey(flatten(i5),
				flatten(i7),
				flatten(p5),
				flatten(p7));
	}
	
	/**
	 * 
	 * @param toFlatten label, potentially with delimiter and index
	 * @return label with delimiter and index removed
	 */
	public static String flatten(String toFlatten){
		if(toFlatten != null && toFlatten.contains(INDEX_SEPARATOR)){
			String [] fields = toFlatten.split(INDEX_SEPARATOR_REGEX);
			return fields[0];
		}
		return toFlatten;
	}
	
	public String getI5Label() {
		return i5;
	}

	public String getI7Label() {
		return i7;
	}

	public String getP5Label() {
		return p5;
	}

	public String getP7Label() {
		return p7;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((i5 == null) ? 0 : i5.hashCode());
		result = prime * result + ((i7 == null) ? 0 : i7.hashCode());
		result = prime * result + ((p5 == null) ? 0 : p5.hashCode());
		result = prime * result + ((p7 == null) ? 0 : p7.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		IndexAndBarcodeKey other = (IndexAndBarcodeKey) obj;
		if (i5 == null) {
			if (other.i5 != null)
				return false;
		} else if (!i5.equals(other.i5))
			return false;
		if (i7 == null) {
			if (other.i7 != null)
				return false;
		} else if (!i7.equals(other.i7))
			return false;
		if (p5 == null) {
			if (other.p5 != null)
				return false;
		} else if (!p5.equals(other.p5))
			return false;
		if (p7 == null) {
			if (other.p7 != null)
				return false;
		} else if (!p7.equals(other.p7))
			return false;
		return true;
	}
	
	
}
