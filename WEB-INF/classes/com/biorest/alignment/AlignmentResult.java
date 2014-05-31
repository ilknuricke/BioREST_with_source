package com.biorest.alignment;

/*
 * Class that encapsulates results of a binary Alignment: parameters, edit distance and the actual alignments
 */
public class AlignmentResult {
	
	private int totalCost=0;
	
	private int alignmentLength=0;
	
	
	public int getAlignmentLength() {
		return alignmentLength;
	}

	public void setAlignmentLength(int alignmentLength) {
		this.alignmentLength = alignmentLength;
	}

	private int SOP=0;
	
	public int getSOP() {
		return SOP;
	}

	public void setSOP(int sOP) {
		SOP = sOP;
	}

	private int matches=0;
	
	private AlignmentParameters parameters=null;
	
	private String[] alignments=null;

	public int getMatches() {
		return matches;
	}

	public void setMatches(int matches) {
		this.matches = matches;
	}

	public int getTotalCost() {
		return totalCost;
	}

	public void setTotalCost(int totalCost) {
		this.totalCost = totalCost;
	}

	public AlignmentParameters getParameters() {
		return parameters;
	}

	public void setParameters(AlignmentParameters parameters) {
		this.parameters = parameters;
	}

	public String[] getAlignments() {
		return alignments;
	}

	public void setAlignments(String[] alignments) {
		this.alignments = alignments;
	}
	
	private int alignmentScore(String seq1, String seq2){
		  totalCost=0;
		  
		  for (int k=0; k < seq1.length(); k++){
			  if (seq1.charAt(k)!=seq2.charAt(k)) {
				   if ( (seq1.charAt(k)!='-') && (seq2.charAt(k)!='-')  ) totalCost++;
			  }			  
			  if ( (seq1.charAt(k)=='-') || (seq2.charAt(k)=='-')  ) {
				  totalCost+=2;
			  }
		  }
		  return totalCost;
	}
	
	public int computeSOP(){
		int _SOP=0;
		for (int i =0; i < alignments.length; i++)
			for (int j =i+1; j < alignments.length; j++){
				_SOP+=alignmentScore(alignments[i],alignments[j]);
			}
		this.SOP=_SOP;
		return SOP;
	}
    
       @Override
        public String toString(){
            
            StringBuilder result = new StringBuilder();
            
            result.append(alignments[0]).append(" size:").append(alignments[0].length()).append("\n");
           
            matches=0;
            int gaps=0;
            for (int k=0; k < alignments[0].length(); k++){
              if (alignments[0].charAt(k)==alignments[1].charAt(k)) {
                   matches++;
                    result.append('|');
             } else result.append(" ");
               
            if ( (alignments[0].charAt(k)=='-') ||
                 (alignments[1].charAt(k)=='-')  )
            gaps++;
            }
           
            result.append("\n");
           
            result.append(alignments[1]).append(" size:").append(alignments[1].length()).append("\n");
            
            float matchRatio=((matches*1.0f)/alignments[0].length())*100;
            System.out.println(matches + " " + matchRatio);
                 
            result.append("Match Score=").append(matches).append(" (").append(matchRatio).append("%)").append(" gaps:").append(gaps).append("\n");
           
             result.append("Edit Distance=").append(getTotalCost()).append("\n");
             result.append("Alignment Length=").append(getAlignmentLength()).append("\n");

             return result.toString();
        }
}
