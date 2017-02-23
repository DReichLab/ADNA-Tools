package adnascreen;

/**
 * Each base-pair has an associated quality. 
 * Quality is stored as a Q-Score [0-40]. 
 * @author Matthew Mah
 *
 */
public class QualitySequence {
	private int[] qScore;
	
	public QualitySequence(String s){
		setQuality(s);
	}
	
	// shallow copy
	public QualitySequence(int[] qScore){
		this.qScore = qScore;
	}
	
	/**
	 * The probability p that a read value is incorrect is
	 * p = 10^(-Q / 10)
	 * @param index
	 * @return
	 */
	public int getQuality(int index){
		return qScore[index];
	}
	
	public void setQuality(String qualityString){
		qScore = new int[qualityString.length()];
		for(int i = 0; i < qualityString.length(); i++){
			qScore[i] = quality33Score(qualityString.charAt(i));
		}
	}
	
	public static char quality33Char(int qScore){
		char c = (char) (qScore + 33);
		return c;
	}
	
	public static int quality33Score(char c){
		return ((int) c) - 33;
	}
	
	public double pIncorrect(int Q){
		double p = Math.pow(10, Q / -10.0); 
		return p;
	}
	
	public int length(){
		return qScore.length;
	}
	
	public QualitySequence subsequence(int beginIndex, int endIndex){
		int size = endIndex - beginIndex;
		if(size < 0)
			throw new IndexOutOfBoundsException();
		int [] scores = new int[size];
		for(int i = 0; i < size; i++){
			scores[i] = qScore[beginIndex + i];
		}
		return new QualitySequence(scores);
	}
	
	@Override
	public String toString(){
		StringBuilder builder = new StringBuilder();
		for(int i = 0; i < qScore.length; i++){
			builder.append(quality33Char(qScore[i]));
		}
		return builder.toString();
	}
	
	public QualitySequence reverse(){
		int [] reversedQScores = new int[qScore.length];
		for (int i = 0; i < qScore.length; i++){
			reversedQScores[i] = qScore[qScore.length - 1 - i];
		}
		return new QualitySequence(reversedQScores);
	}
	
	public boolean equals(Object x){
		if(x instanceof QualitySequence){
			QualitySequence other = (QualitySequence) x;
			if(this.length() == other.length()){
				for(int i = 0; i < this.length(); i++){
					if(this.getQuality(i) != other.getQuality(i)){
						return false;
					}
				}
				return true;
			}
		}
		
		return false;
	}
}
