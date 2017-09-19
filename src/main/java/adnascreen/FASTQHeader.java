package adnascreen;

/**
 * The FASTQ header is the first line of each line quartet. 
	It contains metadata about the instrument and sample parameters.
	This class is for use with Illumina-created FASTQ files.
	Other sources may use different conventions.
 * @author mmah
 *
 */
public class FASTQHeader {
	private String instrument;
	private int runNumber;
	private String flowcellID;
	private int lane;
	private int tile;
	private int x;
	private int y;
	private String UMI;
	private int read;
	private boolean isFiltered;
	private int controlNumber;
	private String index;
	
	// used only for separating fields in getReadGroup results 
	public static final char READ_GROUP_FIELD_DELIMITER = '_';

	public FASTQHeader(String line){
		// htsjdk checks and removes @ header character
		String noLeading; 
		if(line.charAt(0) == '@'){
			noLeading = line.substring(1);
		} else {
			noLeading = line;
		}
		
		String [] initialSplit = noLeading.split(" ");
		if(initialSplit.length > 2)
			throw new IllegalArgumentException("Unexpected spaces in FASTQ header");
		String[] left = initialSplit[0].split(":");
		String[] right = initialSplit[1].split(":");
		
		instrument = left[0];
		runNumber = Integer.valueOf(left[1]);
		flowcellID = left[2];
		lane = Integer.valueOf(left[3]);
		tile = Integer.valueOf(left[4]);
		x = Integer.valueOf(left[5]);
		y = Integer.valueOf(left[6]);
		if(left.length > 7)
			UMI = left[7];
		
		read = Integer.valueOf(right[0]);
		String isFilteredString = right[1];
		switch(isFilteredString){
		case "Y":
			isFiltered = true;
			break;
		case "N":
			isFiltered = false;
			break;
		default:
			throw new IllegalArgumentException("Unexpected value of isFiltered in FASTQ header " + isFilteredString);		
		}
		controlNumber = Integer.valueOf(right[2]);
		index = right[3];
	}
	
	@Override
	public String toString(){
		StringBuilder b = new StringBuilder();
		b.append('@');
		// left side of space
		b.append(instrument);
		b.append(':');
		b.append(runNumber);
		b.append(':');
		b.append(flowcellID);
		b.append(':');
		b.append(lane);
		b.append(':');
		b.append(tile);
		b.append(':');
		b.append(x);
		b.append(':');
		b.append(y);
		if(UMI != null){
			b.append(':');
			b.append(UMI);
		}
		// separating space
		b.append(' ');
		// right side of space
		b.append(read);
		b.append(':');
		b.append(isFiltered ? "Y" : "N");
		b.append(':');
		b.append(controlNumber);
		b.append(':');
		b.append(index);
		
		return b.toString();
	}
	
	/**
	 * This checks whether paired reads and indexed reads belong together based on the FASTQ metadata
	 * @param other
	 * @return
	 */
	public boolean equalsExceptRead(FASTQHeader other){
		return this.instrument.equals(other.instrument)
				&& this.runNumber == other.runNumber
				&& this.flowcellID.equals(other.flowcellID)
				&& this.lane == other.lane
				&& this.tile == other.tile
				&& this.x == other.x
				&& this.y == other.y
				&& (UMI == null || this.UMI.equals(other.UMI))
				// read is not compared
				&& this.isFiltered == other.isFiltered
				&& this.controlNumber == other.controlNumber
				&& this.index.equals(other.index);
	}
	
	@Override
	public boolean equals(Object x){
		if(x instanceof FASTQHeader){
			FASTQHeader other = (FASTQHeader) x;
			return this.equalsExceptRead(other) && (this.read == other.read);
		}
		return false;
	}

	public String getInstrument() {
		return instrument;
	}

	public int getRunNumber() {
		return runNumber;
	}

	public String getFlowcellID() {
		return flowcellID;
	}

	public int getLane() {
		return lane;
	}

	public int getTile() {
		return tile;
	}

	public int getX() {
		return x;
	}

	public int getY() {
		return y;
	}

	public String getUMI() {
		return UMI;
	}

	public int getRead() {
		return read;
	}

	public boolean isFiltered() {
		return isFiltered;
	}

	public int getControlNumber() {
		return controlNumber;
	}

	public String getIndex() {
		return index;
	}
	
	public String getReadGroupElements(){
		StringBuilder b = new StringBuilder();
		b.append("PM:");
		b.append(getInstrument());
		b.append("\t");
		b.append("PU:");
		// platform identifer is flow_cell_id . r . lane
		b.append(getFlowcellID());
		b.append(".");
		b.append(getRunNumber());
		b.append(".");
		b.append(getLane());
		return b.toString();
	}
}
