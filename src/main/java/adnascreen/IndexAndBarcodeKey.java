package adnascreen;

public class IndexAndBarcodeKey {
	private DNASequence i5, i7;
	private String p5BarcodeSet, p7BarcodeSet;
	
	public IndexAndBarcodeKey(DNASequence i5, DNASequence i7, String p5BarcodeSet, String p7BarcodeSet){
		this.i5 = i5;
		this.i7 = i7;
		this.p5BarcodeSet = p5BarcodeSet;
		this.p7BarcodeSet = p7BarcodeSet;
	}
	
	// TODO review this format
	// Consider using JSON, or another more flexible format
	@Override
	public String toString(){
		StringBuilder b = new StringBuilder();
		b.append(i5);
		b.append(',');
		b.append(i7);
		b.append(',');
		b.append(p5BarcodeSet);
		b.append(',');
		b.append(p7BarcodeSet);
		return b.toString();
	}

	public DNASequence getI5() {
		return i5;
	}

	public DNASequence getI7() {
		return i7;
	}

	public String getP5BarcodeSet() {
		return p5BarcodeSet;
	}

	public String getP7BarcodeSet() {
		return p7BarcodeSet;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((i5 == null) ? 0 : i5.hashCode());
		result = prime * result + ((i7 == null) ? 0 : i7.hashCode());
		result = prime * result + ((p5BarcodeSet == null) ? 0 : p5BarcodeSet.hashCode());
		result = prime * result + ((p7BarcodeSet == null) ? 0 : p7BarcodeSet.hashCode());
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
		if (p5BarcodeSet == null) {
			if (other.p5BarcodeSet != null)
				return false;
		} else if (!p5BarcodeSet.equals(other.p5BarcodeSet))
			return false;
		if (p7BarcodeSet == null) {
			if (other.p7BarcodeSet != null)
				return false;
		} else if (!p7BarcodeSet.equals(other.p7BarcodeSet))
			return false;
		return true;
	}
	
	
}
