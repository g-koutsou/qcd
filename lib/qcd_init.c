/* qcd_init.c
 *
 * construction and destruction of structures,
 * definition of geometry,
 *
 * Tomasz Korzec 2008
 **********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

void qcd_antilexic(qcd_uint_2 x[], qcd_uint_8 l, const qcd_uint_2 dim[])
{ 
   l-=(x[0] = (l % dim[0]));   l/=dim[0];
   l-=(x[1] = (l % dim[1]));   l/=dim[1];
   l-=(x[2] = (l % dim[2]));   l/=dim[2];
   x[3]=l;
   return;
}
 
int qcd_initGeometry(qcd_geometry *geo, const qcd_uint_2 L[], const qcd_uint_2 P[], const qcd_real_8 theta[], const int myid, const int numprocs)
{  
   qcd_uint_2 t,x,y,z,b;
   qcd_uint_8 i,j,k;
   qcd_uint_8 lastIndex;

   if( (L[0]%P[0] !=0 ) || (L[1]%P[1] != 0) || (L[2]%P[2] != 0) || (L[3]%P[3] !=0) )
   {
      fprintf(stderr,"process %i: Error in qcd_initGeometry!\n",myid);
      fprintf(stderr,"process %i: %i x %i x %i x %i cannot be subdivided in %i x %i x %i x %i\n ",
              myid,L[0],L[1],L[2],L[3],P[0],P[1],P[2],P[3]);
      return(1);              
   }
   if(P[0]*P[1]*P[2]*P[3] != numprocs)
   {
      fprintf(stderr,"process %i: Error in qcd_initGeometry!\n Need %i processes, got %i\n",myid,P[0]*P[1]*P[2]*P[3],numprocs);
   }
   if(L[0]/P[0] <2 || L[1]/P[1] <2 || L[2]/P[2]<2 || L[3]/P[3]<2)
   {
      fprintf(stderr,"process %i: Error in qcd_initGeometry!\n",myid);
      fprintf(stderr,"process %i: Local lattices must be at least 2x2x2x2\n",myid);
      return(1);
   }
   if( (L[0]/P[0] * L[1]/P[1] * L[2]/P[2]* L[3]/P[3])%2)
   {
      // local volume must be multiple of 2, otherwise E/O partitioning problematic
      fprintf(stderr,"process %i: Error in qcd_initGeometry!\n",myid);
      fprintf(stderr,"process %i: Local volume must be multiple of 2.\n",myid);
      return(1);
   }
   geo->lL[0] = L[0] / P[0];
   geo->lL[1] = L[1] / P[1];
   geo->lL[2] = L[2] / P[2];
   geo->lL[3] = L[3] / P[3];
   geo->lV = geo->lL[0]*geo->lL[1]*geo->lL[2]*geo->lL[3];
   geo->lV3 = geo->lL[1]*geo->lL[2]*geo->lL[3];
   geo->L[0] = L[0];
   geo->L[1] = L[1];
   geo->L[2] = L[2];
   geo->L[3] = L[3];         
   geo->V = geo->L[0]*geo->L[1]*geo->L[2]*geo->L[3];
   geo->V3 = geo->L[1]*geo->L[2]*geo->L[3];
   geo->myid = myid;
   geo->nproc= numprocs;
   geo->theta[0] = theta[0];
   geo->theta[1] = theta[1];
   geo->theta[2] = theta[2];
   geo->theta[3] = theta[3];
    
   qcd_antilexic(geo->Pos, myid, P);
   geo->Pplus[0] = qcd_LEXIC((geo->Pos[0]+1)%P[0],geo->Pos[1],geo->Pos[2],geo->Pos[3],P);
   geo->Pplus[1] = qcd_LEXIC(geo->Pos[0],(geo->Pos[1]+1)%P[1],geo->Pos[2],geo->Pos[3],P);
   geo->Pplus[2] = qcd_LEXIC(geo->Pos[0],geo->Pos[1],(geo->Pos[2]+1)%P[2],geo->Pos[3],P);
   geo->Pplus[3] = qcd_LEXIC(geo->Pos[0],geo->Pos[1],geo->Pos[2],(geo->Pos[3]+1)%P[3],P);
   geo->Pminus[0] = qcd_LEXIC((geo->Pos[0]-1+P[0])%P[0],geo->Pos[1],geo->Pos[2],geo->Pos[3],P);
   geo->Pminus[1] = qcd_LEXIC(geo->Pos[0],(geo->Pos[1]-1+P[1])%P[1],geo->Pos[2],geo->Pos[3],P);
   geo->Pminus[2] = qcd_LEXIC(geo->Pos[0],geo->Pos[1],(geo->Pos[2]-1+P[2])%P[2],geo->Pos[3],P);
   geo->Pminus[3] = qcd_LEXIC(geo->Pos[0],geo->Pos[1],geo->Pos[2],(geo->Pos[3]-1+P[3])%P[3],P);

   geo->plus =  malloc(4*geo->lV*sizeof(qcd_uint_8));
   if(geo->plus == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initGeometry! Out of memory\n",myid);
      return(1);
   }
   geo->minus=  malloc(4*geo->lV*sizeof(qcd_uint_8));
   if(geo->minus == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initGeometry! Out of memory\n",myid);
      return(1);
   }
   
   geo->plus3d =  malloc(4*geo->lV3*sizeof(qcd_uint_8));
   if(geo->plus3d == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initGeometry! Out of memory\n",myid);
      return(1);
   }
   geo->minus3d=  malloc(4*geo->lV3*sizeof(qcd_uint_8));
   if(geo->minus == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initGeometry! Out of memory\n",myid);
      return(1);
   }
   
   for(z=0;z<geo->lL[3];z++)
   for(y=0;y<geo->lL[2];y++)
   for(x=0;x<geo->lL[1];x++)
   for(t=0;t<geo->lL[0];t++)
   {
      i = qcd_LEXIC(t,x,y,z,geo->lL);
      geo->plus[i][0] = qcd_LEXIC((t+1)%geo->lL[0],x,y,z,geo->lL);
      geo->plus[i][1] = qcd_LEXIC(t,(x+1)%geo->lL[1],y,z,geo->lL);
      geo->plus[i][2] = qcd_LEXIC(t,x,(y+1)%geo->lL[2],z,geo->lL);
      geo->plus[i][3] = qcd_LEXIC(t,x,y,(z+1)%geo->lL[3],geo->lL);
      geo->minus[i][0] = qcd_LEXIC((t-1+geo->lL[0])%geo->lL[0],x,y,z,geo->lL);
      geo->minus[i][1] = qcd_LEXIC(t,(x-1+geo->lL[1])%geo->lL[1],y,z,geo->lL);
      geo->minus[i][2] = qcd_LEXIC(t,x,(y-1+geo->lL[2])%geo->lL[2],z,geo->lL);
      geo->minus[i][3] = qcd_LEXIC(t,x,y,(z-1+geo->lL[3])%geo->lL[3],geo->lL);
   }
   lastIndex = geo->lL[0]*geo->lL[1]*geo->lL[2]*geo->lL[3];
   if(geo->lL[0] < geo->L[0]) // are there neighbors in 0 direction?
   {
      for(z=0;z<geo->lL[3];z++)
      for(y=0;y<geo->lL[2];y++)
      for(x=0;x<geo->lL[1];x++)
      {
         i = qcd_LEXIC((geo->lL[0]-1),x,y,z,geo->lL);
         geo->plus[i][0] = qcd_LEXIC0(x,y,z,geo->lL) + lastIndex;
      }
      geo->startBplus[0] = lastIndex;      
      lastIndex += geo->lL[1]*geo->lL[2]*geo->lL[3];
      for(z=0;z<geo->lL[3];z++)
      for(y=0;y<geo->lL[2];y++)
      for(x=0;x<geo->lL[1];x++)
      {
         i = qcd_LEXIC(0,x,y,z,geo->lL);
         geo->minus[i][0] = qcd_LEXIC0(x,y,z,geo->lL) + lastIndex;
      }
      geo->startBminus[0] = lastIndex;
      lastIndex += geo->lL[1]*geo->lL[2]*geo->lL[3];
   }
   if(geo->lL[1] < geo->L[1]) // are there neighbors in 1 direction?
   {
      for(z=0;z<geo->lL[3];z++)
      for(y=0;y<geo->lL[2];y++)
      for(t=0;t<geo->lL[0];t++)
      {
         i = qcd_LEXIC(t,(geo->lL[1]-1),y,z,geo->lL);
         geo->plus[i][1] = qcd_LEXIC1(t,y,z,geo->lL) + lastIndex;
      }
      geo->startBplus[1] = lastIndex;      
      lastIndex += geo->lL[0]*geo->lL[2]*geo->lL[3];
      for(z=0;z<geo->lL[3];z++)
      for(y=0;y<geo->lL[2];y++)
      for(t=0;t<geo->lL[0];t++)
      {
         i = qcd_LEXIC(t,0,y,z,geo->lL);
         geo->minus[i][1] = qcd_LEXIC1(t,y,z,geo->lL) + lastIndex;
      }
      geo->startBminus[1] = lastIndex;
      lastIndex += geo->lL[0]*geo->lL[2]*geo->lL[3];
   }
   if(geo->lL[2] < geo->L[2]) // are there neighbors in 2 direction?
   {
      for(z=0;z<geo->lL[3];z++)
      for(x=0;x<geo->lL[1];x++)
      for(t=0;t<geo->lL[0];t++)
      {
         i = qcd_LEXIC(t,x,(geo->lL[2]-1),z,geo->lL);
         geo->plus[i][2] = qcd_LEXIC2(t,x,z,geo->lL) + lastIndex;
      }
      geo->startBplus[2] = lastIndex;      
      lastIndex += geo->lL[0]*geo->lL[1]*geo->lL[3];
      for(z=0;z<geo->lL[3];z++)
      for(x=0;x<geo->lL[1];x++)
      for(t=0;t<geo->lL[0];t++)
      {
         i = qcd_LEXIC(t,x,0,z,geo->lL);
         geo->minus[i][2] = qcd_LEXIC2(t,x,z,geo->lL) + lastIndex;
      }
      geo->startBminus[2] = lastIndex;
      lastIndex += geo->lL[0]*geo->lL[1]*geo->lL[3];
   }
   if(geo->lL[3] < geo->L[3]) // are there neighbors in 3 direction?
   {
      for(y=0;y<geo->lL[2];y++)
      for(x=0;x<geo->lL[1];x++)
      for(t=0;t<geo->lL[0];t++)
      {
         i = qcd_LEXIC(t,x,y,(geo->lL[3]-1),geo->lL);
         geo->plus[i][3] = qcd_LEXIC3(t,x,y,geo->lL) + lastIndex;
      }
      geo->startBplus[3] = lastIndex;      
      lastIndex += geo->lL[0]*geo->lL[1]*geo->lL[2];
      for(y=0;y<geo->lL[2];y++)
      for(x=0;x<geo->lL[1];x++)
      for(t=0;t<geo->lL[0];t++)
      {
         i = qcd_LEXIC(t,x,y,0,geo->lL);
         geo->minus[i][3] = qcd_LEXIC3(t,x,y,geo->lL) + lastIndex;
      }
      geo->startBminus[3] = lastIndex;
      lastIndex += geo->lL[0]*geo->lL[1]*geo->lL[2];
   } 
   geo->boundarySize=lastIndex-geo->lV;
   
   //tabulate inner boundary points
   geo->edgePoints = geo->lV-(geo->lL[0]-2)*(geo->lL[1]-2)*(geo->lL[2]-2)*(geo->lL[3]-2);
   geo->edge0Points= geo->lV3-(geo->lL[1]-2)*(geo->lL[2]-2)*(geo->lL[3]-2);
   geo->edge = (qcd_uint_8*) malloc(geo->edgePoints*sizeof(qcd_uint_8));
   geo->edge0= (qcd_uint_8*) malloc(geo->edge0Points*sizeof(qcd_uint_8));
   j=0;
   for(z=0; z<geo->lL[3]; z++)
   for(y=0; y<geo->lL[2]; y++)
   for(x=0; x<geo->lL[1]; x++)
   for(t=0; t<geo->lL[0]; t++)
   {
      if(z==0 || z==geo->lL[3]-1 || y==0 || y==geo->lL[2]-1 || x==0 || x==geo->lL[1]-1 || t==0 || t==geo->lL[0]-1)
          geo->edge[j++] = qcd_LEXIC(t,x,y,z,geo->lL);  
   }
   j=0;
   for(z=0; z<geo->lL[3]; z++)
   for(y=0; y<geo->lL[2]; y++)
   for(x=0; x<geo->lL[1]; x++)
   {
      if(z==0 || z==geo->lL[3]-1 || y==0 || y==geo->lL[2]-1 || x==0 || x==geo->lL[1]-1)
          geo->edge0[j++] = qcd_LEXIC0(x,y,z,geo->lL);  
   }
   if(j != geo->edge0Points)
   {
       fprintf(stderr,"Something fishy: %lu inner boundary points instead of %lu\n",j,geo->edge0Points);   
   }


   // neighborhood relations of time-slice objects
   for(z=0;z<geo->lL[3];z++)
   for(y=0;y<geo->lL[2];y++)
   for(x=0;x<geo->lL[1];x++)
   {
      i = qcd_LEXIC0(x,y,z,geo->lL);
      geo->plus3d[i][1] = qcd_LEXIC0((x+1)%geo->lL[1],y,z,geo->lL);
      geo->plus3d[i][2] = qcd_LEXIC0(x,(y+1)%geo->lL[2],z,geo->lL);
      geo->plus3d[i][3] = qcd_LEXIC0(x,y,(z+1)%geo->lL[3],geo->lL);
      geo->minus3d[i][1] = qcd_LEXIC0((x-1+geo->lL[1])%geo->lL[1],y,z,geo->lL);
      geo->minus3d[i][2] = qcd_LEXIC0(x,(y-1+geo->lL[2])%geo->lL[2],z,geo->lL);
      geo->minus3d[i][3] = qcd_LEXIC0(x,y,(z-1+geo->lL[3])%geo->lL[3],geo->lL);
   }
   lastIndex = geo->lL[1]*geo->lL[2]*geo->lL[3];
   if(geo->lL[1] < geo->L[1]) // are there neighbors in 1 direction?
   {
      for(z=0;z<geo->lL[3];z++)
      for(y=0;y<geo->lL[2];y++)
      {
         i = qcd_LEXIC0((geo->lL[1]-1),y,z,geo->lL);
         geo->plus3d[i][1] = qcd_LEXIC01(y,z,geo->lL) + lastIndex;
      }
      geo->startBplus3d[1] = lastIndex;
      lastIndex += geo->lL[2]*geo->lL[3];
      for(z=0;z<geo->lL[3];z++)
      for(y=0;y<geo->lL[2];y++)
      {
         i = qcd_LEXIC0(0,y,z,geo->lL);
         geo->minus3d[i][1] = qcd_LEXIC01(y,z,geo->lL) + lastIndex;
      }
      geo->startBminus3d[1] = lastIndex;
      lastIndex += geo->lL[2]*geo->lL[3];
   }
   if(geo->lL[2] < geo->L[2]) // are there neighbors in 2 direction?
   {
      for(z=0;z<geo->lL[3];z++)
      for(x=0;x<geo->lL[1];x++)
      {
         i = qcd_LEXIC0(x,(geo->lL[2]-1),z,geo->lL);
         geo->plus3d[i][2] = qcd_LEXIC02(x,z,geo->lL) + lastIndex;
      }
      geo->startBplus3d[2] = lastIndex;      
      lastIndex += geo->lL[1]*geo->lL[3];
      for(z=0;z<geo->lL[3];z++)
      for(x=0;x<geo->lL[1];x++)
      {
         i = qcd_LEXIC0(x,0,z,geo->lL);
         geo->minus3d[i][2] = qcd_LEXIC02(x,z,geo->lL) + lastIndex;
      }
      geo->startBminus3d[2] = lastIndex;
      lastIndex += geo->lL[1]*geo->lL[3];
   }
   if(geo->lL[3] < geo->L[3]) // are there neighbors in 3 direction?
   {
      for(y=0;y<geo->lL[2];y++)
      for(x=0;x<geo->lL[1];x++)
      {
         i = qcd_LEXIC0(x,y,(geo->lL[3]-1),geo->lL);
         geo->plus3d[i][3] = qcd_LEXIC03(x,y,geo->lL) + lastIndex;
      }
      geo->startBplus3d[3] = lastIndex;      
      lastIndex += geo->lL[1]*geo->lL[2];
      for(y=0;y<geo->lL[2];y++)
      for(x=0;x<geo->lL[1];x++)
      {
         i = qcd_LEXIC0(x,y,0,geo->lL);
         geo->minus3d[i][3] = qcd_LEXIC03(x,y,geo->lL) + lastIndex;
      }
      geo->startBminus3d[3] = lastIndex;
      lastIndex += geo->lL[1]*geo->lL[2];
   } 
   geo->boundarySize3d=lastIndex-geo->lV3;


   
   // set up partitioning into even/odd sites
   geo->eo=  malloc(geo->lV*sizeof(qcd_uint_8));
   i=0;
   j=0;   
   b=(P[0]*geo->lL[0]+P[1]*geo->lL[1]+P[2]*geo->lL[2]+P[3]*geo->lL[3])%2; // 0: (0,0,0,0) is even, 1: it is odd
   for(z=0; z<geo->lL[3];z++)
   for(y=0; y<geo->lL[2];y++)
   for(x=0; x<geo->lL[1];x++)
   for(t=0; t<geo->lL[0];t++)
   {
      k=qcd_LEXIC(t,x,y,z,geo->lL);
      if( (t+x+y+z+b)%2 )
      {
         geo->eo[i][1]=k;
         i++;
      }else
      {
         geo->eo[j][0]=k;
         j++;
      }
   }
   
   // set up geometry-related MPI stuff
   geo->numOfRequests = 0;
   //4 different derived MPI types for the different vector boundaries
   MPI_Type_vector(geo->lL[1]*geo->lL[2]*geo->lL[3], 4*3*2, geo->lL[0]*4*3*2,MPI_DOUBLE, &(geo->stypeV[0]));
   MPI_Type_vector(geo->lL[2]*geo->lL[3], 4*3*2*geo->lL[0], geo->lL[0]*geo->lL[1]*4*3*2,MPI_DOUBLE, &(geo->stypeV[1]));
   MPI_Type_vector(geo->lL[3], 4*3*2*geo->lL[0]*geo->lL[1], geo->lL[0]*geo->lL[1]*geo->lL[2]*4*3*2,MPI_DOUBLE, &(geo->stypeV[2]));
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*4*3*2, MPI_DOUBLE, &(geo->stypeV[3]));  
   MPI_Type_contiguous(geo->lL[1]*geo->lL[2]*geo->lL[3]*4*3*2, MPI_DOUBLE, &(geo->rtypeV[0]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[2]*geo->lL[3]*4*3*2, MPI_DOUBLE, &(geo->rtypeV[1]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[3]*4*3*2, MPI_DOUBLE, &(geo->rtypeV[2]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*4*3*2, MPI_DOUBLE, &(geo->rtypeV[3]));  
   //and for the different gauge field boundaries
   MPI_Type_vector(geo->lL[1]*geo->lL[2]*geo->lL[3], 4*3*3*2, geo->lL[0]*4*3*3*2,MPI_DOUBLE, &(geo->stypeU[0]));
   MPI_Type_vector(geo->lL[2]*geo->lL[3], 4*3*3*2*geo->lL[0], geo->lL[0]*geo->lL[1]*4*3*3*2,MPI_DOUBLE, &(geo->stypeU[1]));
   MPI_Type_vector(geo->lL[3], 4*3*3*2*geo->lL[0]*geo->lL[1], geo->lL[0]*geo->lL[1]*geo->lL[2]*4*3*3*2,MPI_DOUBLE, &(geo->stypeU[2]));
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*4*3*3*2, MPI_DOUBLE, &(geo->stypeU[3]));  
   MPI_Type_contiguous(geo->lL[1]*geo->lL[2]*geo->lL[3]*4*3*3*2, MPI_DOUBLE, &(geo->rtypeU[0]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[2]*geo->lL[3]*4*3*3*2, MPI_DOUBLE, &(geo->rtypeU[1]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[3]*4*3*3*2, MPI_DOUBLE, &(geo->rtypeU[2]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*4*3*3*2, MPI_DOUBLE, &(geo->rtypeU[3]));  
   //and for the different gauge transformation boundaries
   MPI_Type_vector(geo->lL[1]*geo->lL[2]*geo->lL[3], 3*3*2, geo->lL[0]*3*3*2,MPI_DOUBLE, &(geo->stypeT[0]));
   MPI_Type_vector(geo->lL[2]*geo->lL[3], 3*3*2*geo->lL[0], geo->lL[0]*geo->lL[1]*3*3*2,MPI_DOUBLE, &(geo->stypeT[1]));
   MPI_Type_vector(geo->lL[3], 3*3*2*geo->lL[0]*geo->lL[1], geo->lL[0]*geo->lL[1]*geo->lL[2]*3*3*2,MPI_DOUBLE, &(geo->stypeT[2]));
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*3*3*2, MPI_DOUBLE, &(geo->stypeT[3]));  
   MPI_Type_contiguous(geo->lL[1]*geo->lL[2]*geo->lL[3]*3*3*2, MPI_DOUBLE, &(geo->rtypeT[0]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[2]*geo->lL[3]*3*3*2, MPI_DOUBLE, &(geo->rtypeT[1]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[3]*3*3*2, MPI_DOUBLE, &(geo->rtypeT[2]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*3*3*2, MPI_DOUBLE, &(geo->rtypeT[3]));  
   //and for the different propagator boundaries
   MPI_Type_vector(geo->lL[1]*geo->lL[2]*geo->lL[3], 4*4*3*3*2, geo->lL[0]*4*4*3*3*2,MPI_DOUBLE, &(geo->stypeP[0]));
   MPI_Type_vector(geo->lL[2]*geo->lL[3], 4*4*3*3*2*geo->lL[0], geo->lL[0]*geo->lL[1]*4*4*3*3*2,MPI_DOUBLE, &(geo->stypeP[1]));
   MPI_Type_vector(geo->lL[3], 4*4*3*3*2*geo->lL[0]*geo->lL[1], geo->lL[0]*geo->lL[1]*geo->lL[2]*4*4*3*3*2,MPI_DOUBLE, &(geo->stypeP[2]));
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*4*4*3*3*2, MPI_DOUBLE, &(geo->stypeP[3]));  
   MPI_Type_contiguous(geo->lL[1]*geo->lL[2]*geo->lL[3]*4*4*3*3*2, MPI_DOUBLE, &(geo->rtypeP[0]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[2]*geo->lL[3]*4*4*3*3*2, MPI_DOUBLE, &(geo->rtypeP[1]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[3]*4*4*3*3*2, MPI_DOUBLE, &(geo->rtypeP[2]));  
   MPI_Type_contiguous(geo->lL[0]*geo->lL[1]*geo->lL[2]*4*4*3*3*2, MPI_DOUBLE, &(geo->rtypeP[3]));  
   for(b=0; b<4; b++)
   {
      if(geo->lL[b] < geo->L[b])
      {
         MPI_Type_commit(&(geo->stypeV[b]));
         MPI_Type_commit(&(geo->rtypeV[b]));
         MPI_Type_commit(&(geo->stypeU[b]));
         MPI_Type_commit(&(geo->rtypeU[b]));
         MPI_Type_commit(&(geo->stypeT[b]));
         MPI_Type_commit(&(geo->rtypeT[b]));
         MPI_Type_commit(&(geo->stypeP[b]));
         MPI_Type_commit(&(geo->rtypeP[b]));
      }   
   } 
     
   geo->initialized = 1;
   return(0);
}//end qcd_initGeometry


void qcd_destroyGeometry(const qcd_geometry *geo)
{
   qcd_uint_2 b;
   for(b=0; b<4; b++)
   {
      if(geo->lL[b] < geo->L[b])
      {
         MPI_Type_free((MPI_Datatype*) &(geo->stypeV[b]));
         MPI_Type_free((MPI_Datatype*) &(geo->rtypeV[b]));
         MPI_Type_free((MPI_Datatype*) &(geo->stypeU[b]));
         MPI_Type_free((MPI_Datatype*) &(geo->rtypeU[b]));         
      }
   }   
   free(geo->plus);
   free(geo->minus);
   free(geo->plus3d);
   free(geo->minus3d);
   free(geo->edge);
   free(geo->edge0);
   free(geo->eo);
}//end qcd_destroyGeometry


int qcd_initVector(qcd_vector *vec, qcd_geometry *geo)
{
   qcd_uint_8 size;
   if(!geo->initialized)
   {
      fprintf(stderr,"error in qcd_initVector! Geometry must be properly initialized.\n");
      return(1);
   }
   
   size = geo->lV;
   if(geo->lL[0] < geo->L[0])
      size += 2*(geo->lL[1]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[1] < geo->L[1])
      size += 2*(geo->lL[0]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[2] < geo->L[2])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[3]);
   if(geo->lL[3] < geo->L[3])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[2]);            
   
   vec->D = malloc(4*3*size*sizeof(qcd_complex_16));
   if(vec->D == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initVector! Out of memory\n",geo->myid);
      return(1);
   }
   
   if(geo->lL[0] < geo->L[0])
   {
       vec->Bplus[0]  = &(vec->D[geo->startBplus[0]]  );
       vec->Bminus[0] = &(vec->D[geo->startBminus[0]] );  
   }   
   
   if(geo->lL[1] < geo->L[1])
   {
       vec->Bplus[1]  = &(vec->D[geo->startBplus[1]]);
       vec->Bminus[1] = &(vec->D[geo->startBminus[1]]);  
   }
   
   if(geo->lL[2] < geo->L[2])
   {
       vec->Bplus[2]  = &(vec->D[geo->startBplus[2]]);
       vec->Bminus[2] = &(vec->D[geo->startBminus[2]]);  
   }
   
   if(geo->lL[3] < geo->L[3])
   {
       vec->Bplus[3]  = &(vec->D[geo->startBplus[3]]);
       vec->Bminus[3] = &(vec->D[geo->startBminus[3]]);  
   }
      
   vec->geo = geo;
   vec->initialized = 1;
   return(0);
}//end qcd_initVector

void qcd_destroyVector(qcd_vector *vec)
{
   if(!vec->initialized)
   {
      fprintf(stderr,"Error in qcd_destroyVector! Vector not initialized\n");
      return;
   }
   free(vec->D);
}//end qcd_destroyVector




int qcd_initVector3d(qcd_vector3d *vec, qcd_geometry *geo)
{
   /* initialize structure to sotre a single time-slice of a fermionic field */
   qcd_uint_8 size;
   if(!geo->initialized)
   {
      fprintf(stderr,"error in qcd_initVector3d! Geometry must be properly initialized.\n");
      return(1);
   }
   
   size = geo->lV3;
   if(geo->lL[1] < geo->L[1])
      size += 2*(geo->lL[2]*geo->lL[3]);
   if(geo->lL[2] < geo->L[2])
      size += 2*(geo->lL[1]*geo->lL[3]);
   if(geo->lL[3] < geo->L[3])
      size += 2*(geo->lL[1]*geo->lL[2]);            
   
   vec->D = malloc(4*3*size*sizeof(qcd_complex_16));
   if(vec->D == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initVector3d! Out of memory\n",geo->myid);
      return(1);
   }
      
   if(geo->lL[1] < geo->L[1])
   {
       vec->Bplus[1]  = &(vec->D[geo->startBplus3d[1]]);
       vec->Bminus[1] = &(vec->D[geo->startBminus3d[1]]);  
   }
   
   if(geo->lL[2] < geo->L[2])
   {
       vec->Bplus[2]  = &(vec->D[geo->startBplus3d[2]]);
       vec->Bminus[2] = &(vec->D[geo->startBminus3d[2]]);  
   }
   
   if(geo->lL[3] < geo->L[3])
   {
       vec->Bplus[3]  = &(vec->D[geo->startBplus3d[3]]);
       vec->Bminus[3] = &(vec->D[geo->startBminus3d[3]]);
   }
      
   vec->geo = geo;
   vec->initialized = 1;
   return(0);
}//end qcd_initVector3d

void qcd_destroyVector3d(qcd_vector3d *vec)
{
   if(!vec->initialized)
   {
      fprintf(stderr,"Error in qcd_destroyVector3d! Vector not initialized\n");
      return;
   }
   free(vec->D);
}//end qcd_destroyVector3d



int qcd_initPropagator(qcd_propagator *prp, qcd_geometry *geo)
{
   qcd_uint_8 size;
   if(!geo->initialized)
   {
      fprintf(stderr,"error in qcd_initPropagator! Geometry must be properly initialized.\n");
      return(1);
   }
   
   size = geo->lV;
   if(geo->lL[0] < geo->L[0])
      size += 2*(geo->lL[1]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[1] < geo->L[1])
      size += 2*(geo->lL[0]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[2] < geo->L[2])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[3]);
   if(geo->lL[3] < geo->L[3])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[2]);            

   prp->D = malloc(4*4*3*3*size*sizeof(qcd_complex_16));
   if(prp->D == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initPropagator! Out of memory\n",geo->myid);
      return(1);
   }
   
   if(geo->lL[0] < geo->L[0])
   {
       prp->Bplus[0]  = &(prp->D[geo->startBplus[0]]);
       prp->Bminus[0] = &(prp->D[geo->startBminus[0]]);  
   }   
   
   if(geo->lL[1] < geo->L[1])
   {
       prp->Bplus[1]  = &(prp->D[geo->startBplus[1]]);
       prp->Bminus[1] = &(prp->D[geo->startBminus[1]]);  
   }
   
   if(geo->lL[2] < geo->L[2])
   {
       prp->Bplus[2]  = &(prp->D[geo->startBplus[2]]);
       prp->Bminus[2] = &(prp->D[geo->startBminus[2]]);  
   }
   
   if(geo->lL[3] < geo->L[3])
   {
       prp->Bplus[3]  = &(prp->D[geo->startBplus[3]]);
       prp->Bminus[3] = &(prp->D[geo->startBminus[3]]);  
   }
      
   prp->geo = geo;
   prp->initialized = 1;
   return(0);
}//end qcd_initPropagator

void qcd_destroyPropagator(qcd_propagator *prp)
{
   if(!prp->initialized)
   {
      fprintf(stderr,"Error in qcd_destroyPropagator! Propagator not initialized\n");
      return;
   }
   free(prp->D);
}//end qcd_destroyPropagator




int qcd_initGaugeField(qcd_gaugeField *u, qcd_geometry *geo)
{
   qcd_uint_8 size;
   if(!geo->initialized)
   {
      fprintf(stderr,"error in qcd_initGaugeField! Geometry must be properly initialized.\n");
      return(1);
   }
   size = geo->lV;
   if(geo->lL[0] < geo->L[0])
      size += 2*(geo->lL[1]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[1] < geo->L[1])
      size += 2*(geo->lL[0]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[2] < geo->L[2])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[3]);
   if(geo->lL[3] < geo->L[3])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[2]);            

   u->D = malloc(4*3*3*size*sizeof(qcd_complex_16));
   if(u->D == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initGaugeField! Out of memory\n",geo->myid);
      return(1);
   }
   
   if(geo->lL[0] < geo->L[0])
   {
       u->Bplus[0]  = &(u->D[geo->startBplus[0]]);
       u->Bminus[0] = &(u->D[geo->startBminus[0]]);  
   }   
   
   if(geo->lL[1] < geo->L[1])
   {
       u->Bplus[1]  = &(u->D[geo->startBplus[1]]);
       u->Bminus[1] = &(u->D[geo->startBminus[1]]);  
   }
   
   if(geo->lL[2] < geo->L[2])
   {
       u->Bplus[2]  = &(u->D[geo->startBplus[2]]);
       u->Bminus[2] = &(u->D[geo->startBminus[2]]);  
   }
   
   if(geo->lL[3] < geo->L[3])
   {
       u->Bplus[3]  = &(u->D[geo->startBplus[3]]);
       u->Bminus[3] = &(u->D[geo->startBminus[3]]);  
   }
   u->geo = geo;
   u->initialized = 1;
   return(0);
}//end qcd_initGaugeField


// C. Kallidonis
// Function to set a gauge field to zero
int qcd_setZeroGaugeField(qcd_gaugeField *u, qcd_geometry *geo){

  if(!u->initialized){
    fprintf(stderr,"Error in qcd_setZeroGaugeField! GaugeField not initialized\n");
    return 1;
  }

  for(int v=0;v<geo->lV;v++){
    for(int mu=0;mu<4;mu++){
      for(int c1=0;c1<3;c1++){
	for(int c2=0;c2<3;c2++){
	  u->D[v][mu][c1][c2] = (qcd_complex_16) {0.0,0.0};
	}
      }
    }
  }

  return(0);
}

// C. Kallidonis
// Function to set a gauge field to zero
int qcd_setUnityGaugeField(qcd_gaugeField *u, qcd_geometry *geo){

  if(!u->initialized){
    fprintf(stderr,"Error in qcd_setUnityGaugeField! GaugeField not initialized\n");
    return 1;
  }

  for(int v=0;v<geo->lV;v++){
    for(int mu=0;mu<4;mu++){
      for(int c1=0;c1<3;c1++){
	u->D[v][mu][c1][c1] = (qcd_complex_16) {1.0,0.0};
      }
    }
  }

  return(0);
}


void qcd_destroyGaugeField(qcd_gaugeField *u)
{
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_destroyGaugeField! GaugeField not initialized\n");
      return;
   }
   free(u->D);
}//end qcd_destroyGaugeField


int qcd_initGaugeTransformation(qcd_gaugeTransformation *u, qcd_geometry *geo)
{
   qcd_uint_8 size;
   if(!geo->initialized)
   {
      fprintf(stderr,"error in qcd_initGaugeField! Geometry must be properly initialized.\n");
      return(1);
   }
   size = geo->lV;
   if(geo->lL[0] < geo->L[0])
      size += 2*(geo->lL[1]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[1] < geo->L[1])
      size += 2*(geo->lL[0]*geo->lL[2]*geo->lL[3]);
   if(geo->lL[2] < geo->L[2])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[3]);
   if(geo->lL[3] < geo->L[3])
      size += 2*(geo->lL[0]*geo->lL[1]*geo->lL[2]);            

   u->D = malloc(3*3*size*sizeof(qcd_complex_16));
   if(u->D == NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_initGaugeTransformation! Out of memory\n",geo->myid);
      return(1);
   }
   
   if(geo->lL[0] < geo->L[0])
   {
       u->Bplus[0]  = &(u->D[geo->startBplus[0]]);
       u->Bminus[0] = &(u->D[geo->startBminus[0]]);  
   }   
   
   if(geo->lL[1] < geo->L[1])
   {
       u->Bplus[1]  = &(u->D[geo->startBplus[1]]);
       u->Bminus[1] = &(u->D[geo->startBminus[1]]);  
   }
   
   if(geo->lL[2] < geo->L[2])
   {
       u->Bplus[2]  = &(u->D[geo->startBplus[2]]);
       u->Bminus[2] = &(u->D[geo->startBminus[2]]);  
   }
   
   if(geo->lL[3] < geo->L[3])
   {
       u->Bplus[3]  = &(u->D[geo->startBplus[3]]);
       u->Bminus[3] = &(u->D[geo->startBminus[3]]);  
   }
   u->geo = geo;
   u->initialized = 1;
   return(0);
}//end qcd_initGaugeTransformation


void qcd_destroyGaugeTransformation(qcd_gaugeTransformation *u)
{
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_destroyGaugeTransformation! GaugeTransformation not initialized\n");
      return;
   }
   free(u->D);
}//end qcd_destroyGaugeTransformation
