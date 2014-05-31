package com.biorest.alignment;

public class NeedlemanWunsch3D {
	public static boolean verbose=false;
	
	
	/*
	 * prints out the score matrix 
	 */
	public static void dumpMatrix(int[][] matrix, String row, String column){
	    
		 System.out.print(String.format("%5s",""));
		 for (int j =0; j< row.length(); j++){
		   System.out.print(String.format("%5s", row.charAt(j)+"  "));
		 }	    	
		 System.out.println();
		 
		 for (int i =0; i< column.length(); i++){
			   System.out.print(String.format("%5s",column.charAt(i)+ " "));
		    	for (int j =0; j< row.length(); j++){
		    		System.out.print(String.format("%5s", matrix[i][j]+" "));
		    	}
		    	System.out.println();
		   }
	 }
	
	 public static int[][] NeedlemanWunsch2D(String seq1, String seq2, SimpleAlignmentParameters parameters) {
			
	int[][] scoreMatrix = new int[seq1.length()+1][seq2.length()+1];
    
    int gapPenalty=parameters.getGapPenalty();
    int substitutePenalty=parameters.getSubstitutePenalty();
    
    //Initialize the score matrix
    //the first row and column are for the gap
    //Complexity: O(NxM)
    for (int i =0; i< seq1.length()+1; i++){
    	for (int j =0; j< seq2.length()+1; j++){
	     scoreMatrix[i][j]=0;
	     if (i==0) 
	    	 scoreMatrix[i][j] = gapPenalty*j;
	     else if (j==0) 
	    	 scoreMatrix[i][j] = gapPenalty*i;
	 }
    }


    int similarityCost=0;
    int matchCost=0;
    int seq1GapCost=0;
    int seq2GapCost=0;
    
    //Compute the minimum cost scores between all possible pairs of prefixes
    //Complexity: O(NxM)
    for (int i =1; i< seq1.length()+1; i++){
    	for (int j =1; j< seq2.length()+1; j++){
    		
    		//Case 1: The cost of mistmatch between the two prefixes
    		similarityCost= (seq2.charAt(j-1)==seq1.charAt(i-1)) ? 0 : substitutePenalty;   
    		matchCost = scoreMatrix[i-1][j-1] + similarityCost;
    		
    		//Case 2: the cost of adding a gap on sequence 2
    		seq2GapCost = scoreMatrix[i-1][j] + gapPenalty;
    		
    		//Case 3: the cost of adding a gap on sequence 1
    		seq1GapCost = scoreMatrix[i][j-1] + gapPenalty;
    		
    		//System.out.println(matchCost + " " + deleteCost + " " + insertionCost);
    		
    		scoreMatrix[i][j] = Math.min(Math.min(matchCost,seq1GapCost),seq2GapCost);
    	}
    }
    
     return scoreMatrix;
   }
	 
	 /*
	  * Needleman-Wunsch Dynamic Programming Algorithm
	  * see http://amrita.vlab.co.in/?sub=3&brch=274&sim=1431&cnt=1
	  * Runtime complexity: O(NxMxK), Space complexity: O(NxMxK) where N,M,K are lengths of the sequences
	  */
	 public static AlignmentResult computeNWAlignment(String seq1, String seq2, String seq3,SimpleAlignmentParameters parameters) {
			

		 int[][][] scoreMatrix = new int[seq1.length()+1][seq2.length()+1][seq3.length()+1];
		    
		    int gapPenalty=parameters.getGapPenalty();
		    int substitutePenalty=parameters.getSubstitutePenalty();
		    
		    AlignmentResult result = new AlignmentResult();
		    result.setParameters(parameters);
		    
		    int[][] ijscoreMatrix = NeedlemanWunsch2D(seq1, seq2, parameters);
		    int[][] jkscoreMatrix = NeedlemanWunsch2D(seq2, seq3, parameters);
		    int[][] ikscoreMatrix = NeedlemanWunsch2D(seq1, seq3, parameters);
		    
		    //Initialize the score matrix
		    //the first rows of each dimension
		    //Complexity: O(NxMxK)
		    for (int i =0; i< seq1.length()+1; i++){
			 for (int j =0; j< seq2.length()+1; j++){
		    	for (int k =0; k< seq3.length()+1; k++){
			     scoreMatrix[i][j][k]=0;
			     //1D
			     if (i==0&&k==0) 
			    	 scoreMatrix[i][j][k] = gapPenalty*j*2;
			     else if (i==0&&j==0) 
			    	 scoreMatrix[i][j][k] = gapPenalty*k*2;
			     else if (j==0&&k==0) 
			    	 scoreMatrix[i][j][k] = gapPenalty*i*2;
			    //2D
			     else if (i==0)
			     	 scoreMatrix[i][j][k] = jkscoreMatrix[j][k]+ gapPenalty*(j+k);
			     else if (j==0)
			     	 scoreMatrix[i][j][k] = ikscoreMatrix[i][k]+ gapPenalty*(i+k);
			     else if (k==0)
			     	 scoreMatrix[i][j][k] = ijscoreMatrix[i][j]+ gapPenalty*(i+j);
			     
			     System.out.print(i + " " + j + " " +k + " : "+scoreMatrix[i][j][k]+"\n " );
		    	}
		     }
		    }
		    
		    NeedlemanWunsch.dumpMatrix(ijscoreMatrix, " " + seq1 , " " + seq2);
		    NeedlemanWunsch.dumpMatrix(ikscoreMatrix, " " + seq1 , " " + seq3);
		    NeedlemanWunsch.dumpMatrix(jkscoreMatrix, " " + seq2 , " " + seq3);
			    
		    int similarityCost=0;
		    int matchCost=0;
		    int seq1GapCost=0;
		    int seq2GapCost=0;
		    int seq3GapCost=0;
		    
		    int seq12Cost=0;
		    int seq13Cost=0;
		    int seq23Cost=0;
		    
		    //Compute the minimum cost scores between all possible pairs of prefixes
		    //Complexity: O(NxMxK)
		    for (int i =1; i< seq1.length()+1; i++){
				 for (int j =1; j< seq2.length()+1; j++){
			    	for (int k =1; k< seq3.length()+1; k++){
						
			    	//Now there will be 7 cases
			    		
		    		//Case 1: The cost of mistmatch between the three prefixes
			        similarityCost=0;
		    		similarityCost+= (seq3.charAt(k-1)==seq2.charAt(j-1)) ? 0 : substitutePenalty;   
		    		similarityCost+= (seq3.charAt(k-1)==seq1.charAt(i-1)) ? 0 : substitutePenalty;   
		    		similarityCost+= (seq2.charAt(j-1)==seq1.charAt(i-1)) ? 0 : substitutePenalty;   
		    		matchCost = scoreMatrix[i-1][j-1][k-1] + similarityCost;
		    		
		    		//Case 2: 
		    		seq1GapCost = scoreMatrix[i-1][j][k] + gapPenalty*2;
		    		
		    		//Case 3: 
		    		seq2GapCost = scoreMatrix[i][j-1][k] + gapPenalty*2;
		    		
		    		//Case 4: 
		    		seq3GapCost = scoreMatrix[i][j][k-1] + gapPenalty*2;
		    		
		    		//Case 5 : 
		    		seq12Cost  = scoreMatrix[i-1][j-1][k] + gapPenalty + ((seq1.charAt(i-1)==seq2.charAt(j-1)) ? 0 : substitutePenalty);
		    			
		    		//Case 6 : 
		    		seq13Cost  = scoreMatrix[i-1][j][k-1] + gapPenalty + ((seq1.charAt(i-1)==seq3.charAt(k-1)) ? 0 : substitutePenalty);
		    		//System.out.println("** " +scoreMatrix[i-1][j][k-1] + "   " + seq1.charAt(i-1) + "   " + seq3.charAt(k-1) + "  "+ ((seq1.charAt(i-1)==seq3.charAt(k-1)) ? 0 : substitutePenalty)  );
		    		
		    		//Case 7 :  
		    		seq23Cost  = scoreMatrix[i][j-1][k-1] + gapPenalty + ((seq2.charAt(j-1)==seq3.charAt(k-1)) ? 0 : substitutePenalty);
		    		
		    		System.out.println(matchCost + " " + seq1GapCost + " " + seq2GapCost + " " + seq3GapCost+ " " + seq12Cost + " " +seq13Cost+" "  + seq23Cost);
		    		
		    		scoreMatrix[i][j][k] = Math.min(Math.min(matchCost,seq1GapCost),
		    				                        Math.min(Math.min(Math.min(seq2GapCost,seq3GapCost),
		    				                        Math.min(seq13Cost, seq23Cost)),seq12Cost));
		    		  System.out.print(scoreMatrix[i][j][k]+ " \n");
			    	}
			  }
		    }
		    
		    //Reconstruct the Alignment by backtracking on the score matrix
		    //Complexity O(N)
		    StringBuilder alignedSequence1= new StringBuilder();
		    StringBuilder alignedSequence2= new StringBuilder();
		    StringBuilder alignedSequence3= new StringBuilder();
		    
		    int i = seq1.length();
		    int j = seq2.length();
		    int k = seq3.length();
		    
		    while (i >0 || j > 0 || k > 0) {
		    
		       	System.out.println(i+" " + j + " " + k + ": " + scoreMatrix[i][j][k]);
				 
		    	if (i >0 && j > 0 && k>0) {
		    		similarityCost=0;
		    		similarityCost+= (seq3.charAt(k-1)==seq2.charAt(j-1)) ? 0 : substitutePenalty;   
		    		similarityCost+= (seq3.charAt(k-1)==seq1.charAt(i-1)) ? 0 : substitutePenalty;   
		    		similarityCost+= (seq2.charAt(j-1)==seq1.charAt(i-1)) ? 0 : substitutePenalty;   
		    	}
		    	else similarityCost = Integer.MAX_VALUE;
		    	
		    	//Case 1
		    	if ( i > 0 && j > 0 && k > 0 && (scoreMatrix[i][j][k] == scoreMatrix[i-1][j-1][k-1] + similarityCost)) { 
		    		alignedSequence1.append(seq1.charAt(i-1));
		    		alignedSequence2.append(seq2.charAt(j-1));
		    		alignedSequence3.append(seq3.charAt(k-1));
		    		i=i-1;
		    		j=j-1;
		    		k=k-1;
		    		System.out.println(" -->" + 1);
		    	}
		    	else if ( i> 0 && j>0 && scoreMatrix[i][j][k] == (scoreMatrix[i-1][j-1][k] + gapPenalty +((seq1.charAt(i-1)==seq2.charAt(j-1)) ? 0 : substitutePenalty))){ //case 5
		    		alignedSequence3.append("-");
		    		alignedSequence2.append(seq2.charAt(j-1));
		    		alignedSequence1.append(seq1.charAt(i-1));
		    		i=i-1;
		    		j=j-1;
		    		System.out.println(" -->" + 2);
			    	
		    	}else if ( i > 0 && k > 0 && scoreMatrix[i][j][k] == (scoreMatrix[i-1][j][k-1] + gapPenalty+((seq3.charAt(k-1)==seq1.charAt(i-1)) ? 0 : substitutePenalty))){ //case 6
		    		alignedSequence3.append(seq3.charAt(k-1));
		    		alignedSequence2.append("-");
		    		alignedSequence1.append(seq1.charAt(i-1));
		    		i=i-1;
		    		k=k-1;
		    		System.out.println(" -->" + 3);
			    	
		    	}else if ( j > 0 && k > 0 && scoreMatrix[i][j][k] == (scoreMatrix[i][j-1][k-1] + gapPenalty+((seq2.charAt(j-1)==seq3.charAt(k-1)) ? 0 : substitutePenalty))){ //case 7
		    		alignedSequence3.append(seq3.charAt(k-1));
		    		alignedSequence2.append(seq2.charAt(j-1));
		    		alignedSequence1.append("-");
		    		j=j-1;
		    		k=k-1;
		    		System.out.println(" -->" + 4);
			    	
		    	}else if ( i > 0 && scoreMatrix[i][j][k] == (scoreMatrix[i-1][j][k] + gapPenalty*2)){ //case 2
		    		alignedSequence3.append("-");
		    		alignedSequence2.append("-");
		    		alignedSequence1.append(seq1.charAt(i-1));
		    		i=i-1;
		    		System.out.println(" -->" + 5);
			    	
		    	}
		    	else if ( j > 0 && scoreMatrix[i][j][k] == (scoreMatrix[i][j-1][k] + gapPenalty*2)){ //case 3
		    		alignedSequence3.append("-");
		    		alignedSequence2.append(seq2.charAt(j-1));
		    		alignedSequence1.append("-");
		    		j=j-1;
		    		System.out.println(" -->" + 6);
			    	
		    	}
		    	else if ( k > 0 && scoreMatrix[i][j][k] == (scoreMatrix[i][j][k-1] + gapPenalty*2)){ //case 4
		    		alignedSequence3.append(seq3.charAt(k-1));
		    		alignedSequence2.append("-");
		    		alignedSequence1.append("-");
		    		k=k-1;
		    		System.out.println(" -->" + 7);
			    	
		    	}else{
		    		return null; 
		    	}
		    	
		    } // end of while
		    
		    
		    result.setTotalCost(scoreMatrix[seq1.length()][seq2.length()][seq3.length()]);
		    result.setAlignmentLength(alignedSequence1.length());
		    result.setAlignments(new String[] {alignedSequence1.reverse().toString(), 
		    								   alignedSequence2.reverse().toString(),
		    								   alignedSequence3.reverse().toString()});
		 	    
		   return result;
		     }
}

 