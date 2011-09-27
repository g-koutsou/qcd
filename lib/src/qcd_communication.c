/* qcd_communication.c
 * 
 * some non-blocking communication routines
 * use: MPI
 *
 * Tomasz Korzec 2008
 **********************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
 
/* start sending vector-field boundaries to + and - neighbors
 * and gauge-field boundaries to + neighbors
 * that's what is needed in wilson dirac operators
 */ 
int qcd_communicateVectorPMGaugeP(qcd_vector *v, qcd_gaugeField *u)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(v->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateVectorPMGaugeP! Previous communication not finished\n",v->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(v->geo->lL[b] < v->geo->L[b])
      {
         //start communication in b direction
         //send to - / recieve from +
         MPI_Isend(&(v->D[0][0][0]), 1, v->geo->stypeV[b], v->geo->Pminus[b], 3*b, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(v->Bplus[b][0][0]), 1, v->geo->rtypeV[b], v->geo->Pplus[b], 3*b, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         
         startpos=v->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= v->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(v->D[startpos][0][0]), 1, v->geo->stypeV[b], v->geo->Pplus[b], 3*b+1, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(v->Bminus[b][0][0]), 1, v->geo->rtypeV[b], v->geo->Pminus[b], 3*b+1, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Isend(&(u->D[startpos][0][0][0]), 1, v->geo->stypeU[b], v->geo->Pplus[b], 3*b+2, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(u->Bminus[b][0][0][0]), 1, v->geo->rtypeU[b], v->geo->Pminus[b], 3*b+2, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicateVectorPMGaugeP


/* start sending vector-field boundaries to + and - neighbors
 */ 
int qcd_communicateVectorPM(qcd_vector *v)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(v->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateVectorPM! Previous communication not finished\n",v->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(v->geo->lL[b] < v->geo->L[b])
      {
         //start communication in b direction
         //send to - / recieve from +
         MPI_Isend(&(v->D[0][0][0]), 1, v->geo->stypeV[b], v->geo->Pminus[b], 3*b, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(v->Bplus[b][0][0]), 1, v->geo->rtypeV[b], v->geo->Pplus[b], 3*b, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         
         startpos=v->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= v->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(v->D[startpos][0][0]), 1, v->geo->stypeV[b], v->geo->Pplus[b], 3*b+1, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
         MPI_Irecv(&(v->Bminus[b][0][0]), 1, v->geo->rtypeV[b], v->geo->Pminus[b], 3*b+1, MPI_COMM_WORLD, &(v->geo->requests[v->geo->numOfRequests++]));
      }
      
   }//end communication-start loop
      
   return 0;
}//end qcd_communicateVectorPM



/* start sending gauge-field boundaries to - neighbors
 * that's what is needed to calculate the plaquette
 */ 
int qcd_communicateGaugeM(qcd_gaugeField *u)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(u->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateGaugeM! Previous communication not finished\n",u->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(u->geo->lL[b] < u->geo->L[b])
      {
         //start communication in b direction
         //send to - / recieve from +
         MPI_Isend(&(u->D[0][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pminus[b], 3*b, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
         MPI_Irecv(&(u->Bplus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pplus[b], 3*b, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));         
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicateGaugeM


/* start gauge-field boundaries to + neighbors
 */ 
int qcd_communicateGaugeP(qcd_gaugeField *u)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(u->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateGaugeP! Previous communication not finished\n",u->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(u->geo->lL[b] < u->geo->L[b])
      {         
         startpos=u->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= u->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(u->D[startpos][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pplus[b], 3*b+2, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
         MPI_Irecv(&(u->Bminus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pminus[b], 3*b+2, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicateGaugeP



/* "Inverse" communication: send _outer_ boundaries to + neighbor's inner boundary 
 */ 
int qcd_communicateGaugePInverse(qcd_gaugeField *u)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(u->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateGaugePInverse! Previous communication not finished\n",u->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(u->geo->lL[b] < u->geo->L[b])
      {         
         startpos=u->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= u->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(u->Bminus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pplus[b], 3*b+2, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
         MPI_Irecv(&(u->D[startpos][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pminus[b], 3*b+2, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicateGaugePInverse




/* start sending gauge-field boundaries to + and - neighbors
 * that's what is needed during APE-smearing
 */ 
int qcd_communicateGaugePM(qcd_gaugeField *u)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(u->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateGaugePM! Previous communication not finished\n",u->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(u->geo->lL[b] < u->geo->L[b])
      {
         //start communication in b direction
         //send to - / recieve from +
         MPI_Isend(&(u->D[0][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pminus[b], 3*b, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
         MPI_Irecv(&(u->Bplus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pplus[b], 3*b, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
         
         startpos=u->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= u->geo->lL[bb];
            
         //send to + / recieve from -
         //send to + / recieve from -
         MPI_Isend(&(u->D[startpos][0][0][0]), 1, u->geo->stypeU[b], u->geo->Pplus[b], 3*b+2, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
         MPI_Irecv(&(u->Bminus[b][0][0][0]), 1, u->geo->rtypeU[b], u->geo->Pminus[b], 3*b+2, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicateGaugePM


/* start sending gauge-transformation field boundaries to - neighbors
 * that's what is needed to  carry out a gauge transformation
 */ 
int qcd_communicateTransformationM(qcd_gaugeTransformation *u)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(u->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicateTransformationM! Previous communication not finished\n",u->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(u->geo->lL[b] < u->geo->L[b])
      {
         //start communication in b direction
         //send to - / recieve from +
         MPI_Isend(&(u->D[0][0][0]), 1, u->geo->stypeT[b], u->geo->Pminus[b], 3*b, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
         MPI_Irecv(&(u->Bplus[b][0][0]), 1, u->geo->rtypeT[b], u->geo->Pplus[b], 3*b, MPI_COMM_WORLD, &(u->geo->requests[u->geo->numOfRequests++]));
      }      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicateTransformationM



/* start sending propagator boundaries to + neighbors
 * needed for certain observables
 */ 
int qcd_communicatePropagatorP(qcd_propagator *p)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(p->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicatePropagatorP! Previous communication not finished\n",p->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(p->geo->lL[b] < p->geo->L[b])
      {
      
         startpos=p->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= p->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(p->D[startpos][0][0][0][0]), 1, p->geo->stypeP[b], p->geo->Pplus[b], 3*b+1, MPI_COMM_WORLD, &(p->geo->requests[p->geo->numOfRequests++]));
         MPI_Irecv(&(p->Bminus[b][0][0][0][0]), 1, p->geo->rtypeP[b], p->geo->Pminus[b], 3*b+1, MPI_COMM_WORLD, &(p->geo->requests[p->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicatePropagatorP



int qcd_communicatePropagatorPM(qcd_propagator *p)
{
   qcd_uint_2 b,bb;
   qcd_uint_8 startpos;
   
   if(p->geo->numOfRequests != 0)
   {
      fprintf(stderr,"process %i: Error in qcd_communicatePropagatorPM! Previous communication not finished\n",p->geo->myid);
      return(1);
   }
   
   //start communication   
   for(b=0; b<4; b++)
   {
      if(p->geo->lL[b] < p->geo->L[b])
      {
         //start communication in b direction
         //send to - / recieve from +
         MPI_Isend(&(p->D[0][0][0][0][0]), 1, p->geo->stypeP[b], p->geo->Pminus[b], 3*b, MPI_COMM_WORLD, &(p->geo->requests[p->geo->numOfRequests++]));
         MPI_Irecv(&(p->Bplus[b][0][0][0][0]), 1, p->geo->rtypeP[b], p->geo->Pplus[b], 3*b, MPI_COMM_WORLD, &(p->geo->requests[p->geo->numOfRequests++]));
         
         startpos=p->geo->lL[b]-1;
         for(bb=0; bb<b; bb++)
            startpos *= p->geo->lL[bb];
            
         //send to + / recieve from -
         MPI_Isend(&(p->D[startpos][0][0][0][0]), 1, p->geo->stypeP[b], p->geo->Pplus[b], 3*b+1, MPI_COMM_WORLD, &(p->geo->requests[p->geo->numOfRequests++]));
         MPI_Irecv(&(p->Bminus[b][0][0][0][0]), 1, p->geo->rtypeP[b], p->geo->Pminus[b], 3*b+1, MPI_COMM_WORLD, &(p->geo->requests[p->geo->numOfRequests++]));
      }
      
   }//end communication-start loop  
      
   return 0;
}//end qcd_communicatePropagatorPM





void qcd_waitall(qcd_geometry *geo)
{
   //if(geo->myid==0) printf("waiting for %i requests to finish\n",geo->numOfRequests);
   MPI_Waitall(geo->numOfRequests, geo->requests, MPI_STATUSES_IGNORE);
   geo->numOfRequests=0;
}//end qcd_waitall
