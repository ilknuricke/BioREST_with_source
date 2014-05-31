
package com.biorest.rest;

import javax.ejb.EJB;
import javax.ejb.Stateless;

import javax.ws.rs.*;
import javax.ws.rs.core.Response;




@Stateless
@Path("/alignment")
public class AlignmentWS {

    @EJB
    private AlignmentBean alignmentBean;
    @GET
    @Produces("text/plain")
    public String getXml(@DefaultValue("") @QueryParam("sequences") String sequences) {
             System.out.println("sequences= "+sequences);
            
             return alignmentBean.getAlignment(sequences);
    }

    /**
     * PUT method for updating an instance of HelloWorldResource
     * @param content representation for the resource
     * @return an HTTP response with content of the updated or created resource.
     */
    @PUT
    @Consumes("text/plain")
    public void putXml(String content) {
        //nothing
    }
    
     
   @OPTIONS
   public Response getOptions() {
       
    return Response.ok()
      .header("Access-Control-Allow-Origin", "*")
      .header("Access-Control-Allow-Methods", "POST, GET, PUT, UPDATE, OPTIONS")
      .header("Access-Control-Allow-Headers", "Content-Type, Accept, X-Requested-With").build();
  }
}
