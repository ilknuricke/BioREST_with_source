
package com.biorest.rest;

import com.biorest.alignment.AlignmentResult;
import com.biorest.alignment.NeedlemanWunsch;
import com.biorest.alignment.SimpleAlignmentParameters;
import javax.ejb.Singleton;

/**
 *
 * @author ilknur
 */
@Singleton
public class AlignmentBean {
    public String getAlignment (String sequences){
        
       String[] seqs=sequences.split(",");
       AlignmentResult result = NeedlemanWunsch.computeNWAlignment(seqs[0], seqs[1], new SimpleAlignmentParameters());
        
       
      return result.toString();
    }
}
