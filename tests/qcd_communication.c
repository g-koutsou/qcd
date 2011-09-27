/* qcd_communication.c
 * 
 * some non-blocking communication routines
 * use: MPI
 *
 * Tomasz Korzec 2008
 **********************************************/
 
 
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


void qcd_waitall(qcd_geometry *geo)
{
   //if(geo->myid==0) printf("waiting for %i requests to finish\n",geo->numOfRequests);
   MPI_Waitall(geo->numOfRequests, geo->requests, MPI_STATUSES_IGNORE);
   geo->numOfRequests=0;
}//end qcd_waitall
