package com.biorest.alignment;

/*
 *  Alignment parameters used in Needleman-Wunsch Algorithm with gap penalty
 */
public class SimpleAlignmentParameters extends AlignmentParameters{

	private int gapPenalty;
	
	private int substitutePenalty;

	public int getGapPenalty() {
		return gapPenalty;
	}

	public void setGapPenalty(int gapPenalty) {
		this.gapPenalty = gapPenalty;
	}

	public int getSubstitutePenalty() {
		return substitutePenalty;
	}

	public void setSubstitutePenalty(int substitutePenalty) {
		this.substitutePenalty = substitutePenalty;
	}

	public SimpleAlignmentParameters(int gapPenalty, int substitutePenalty) {
		super();
		this.gapPenalty = gapPenalty;
		this.substitutePenalty = substitutePenalty;
	}
	
	public SimpleAlignmentParameters() {
		super();
		this.gapPenalty = 2;
		this.substitutePenalty = 1;
	}
}
