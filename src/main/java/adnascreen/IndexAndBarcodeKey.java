package adnascreen;

public class IndexAndBarcodeKey {
	// labels
	private String i5, i7;
	private String p5, p7;

	// - is chosen because it is a natural character in filenames 
	public static final char FIELD_SEPARATOR = '-'; // this needs to be distinct from MergedRead.KEY_SEPARATOR
	
	public IndexAndBarcodeKey(String i5, String i7, String p5, String p7){
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
		p5 = fields[2];
		p7 = fields[3];
	}
	
	// TODO review this format
	// Consider using JSON, or another more flexible format
	@Override
	public String toString(){
		StringBuilder b = new StringBuilder();
		b.append(i5);
		b.append(FIELD_SEPARATOR);
		b.append(i7);
		b.append(FIELD_SEPARATOR);
		b.append(p5);
		b.append(FIELD_SEPARATOR);
		b.append(p7);
		return b.toString();
	}

	public String getI5() {
		return i5;
	}

	public String getI7() {
		return i7;
	}

	public String getP5() {
		return p5;
	}

	public String getP7() {
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
