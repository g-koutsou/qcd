/* qcd_io.c
 *
 * parallel I/O routines
 * for gauge fields and 
 * propagators in several
 * formats.
 *
 * Tomasz Korzec 2008
 ******************************************/

#include <lime.h>

/*
  Function that returns 1 if machine is big endian,
  0 if little endina. (by Giannis Koutsou)
*/ 
int qcd_isBigEndian()
{
   union{
     char C[4];
     int  R   ;
        }word;
   word.R=1;
   if(word.C[3]==1) return 1;
   if(word.C[0]==1) return 0;
}


/*
  Function that converts between
  big endian and little endian
  8 byte words.
*/
void qcd_swap_8(double *Rd, int N)
{
   register char *i,*j,*k;
   char swap;
   char *max;
   char *R = (char*) Rd;

   max = R+(N<<3);
   for(i=R;i<max;i+=8)
   {
      j=i; k=j+7;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
      j++; k--;
      swap = *j; *j = *k;  *k = swap;
   }
}


void qcd_swap_4(float *Rd, int N)
{
  register char *i,*j,*k;
  char swap;
  char *max;
  char *R =(char*) Rd;

  max = R+(N<<2);
  for(i=R;i<max;i+=4)
  {
    j=i; k=j+3;
    swap = *j; *j = *k;  *k = swap;
    j++; k--;
    swap = *j; *j = *k;  *k = swap;
  }
}


/* by Giannis Koutsou */
char* qcd_getParams(char* fname,int *len)
{
   FILE *pfile;
   char *params;
   int i;
  
   if ((pfile=fopen(fname,"r"))==NULL)
   {
      fprintf(stderr,"Error, cannot open %s for reading\n",fname);
      return(NULL);
   }
  
   i=0;
   while(!feof(pfile))
   {
      fread(&params,1,1,pfile);
      i++;
   }
   *(len)=i;
   rewind(pfile);
   params = (char*)malloc(*(len)*sizeof(char));

   for(i=0;i<*(len);i++)
      fread(params+i,1,1,pfile);
   
   fclose(pfile);
  
   return params;
}


/* by Giannis Koutsou */
char* qcd_getParam(char token[],char* params,int len)
{
   int i,token_len=strlen(token);
   char*value;

   for(i=0;i<len;i++)
   {
      if(memcmp(token,params+i,token_len)==0)
      {
         i+=token_len;
         *(index(params+i,'<'))='\0';
         break;
      }
   }
   return params+i;
}
 

 
int qcd_getGaugeLime(char *fname, qcd_gaugeField *u)
{
   LimeReader *limereader;
   FILE *fid;
   char *lime_type,*lime_data;
   n_uint64_t lime_data_size;
   char dummy;
   MPI_Offset offset;
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   int sizes[5], lsizes[5], starts[5];
   qcd_uint_8 i,j;
   qcd_uint_4 stride;
   qcd_uint_2 chunksize,mu,nu,c1,c2;
   char *swap[16*3*3];
   char *buffer;
   qcd_uint_4 x,y,z,t;
   int  isDouble;
      
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_getGaugeLime! Gauge field must be properly initialized.\n");
      return(1);
   }   
   
   if(u->geo->myid == 0)
   {
      /* read lime header */
      fid=fopen(fname,"r");
      if(fid==NULL)
      {
         fprintf(stderr,"process 0: Error in qcd_getGaugeLime! Could not open %s for reading\n",fname);
         return(1);
      }
      if ((limereader = limeCreateReader(fid))==NULL)
      {
         fprintf(stderr,"process 0: Error in qcd_getGaugeLime! Could not create limeReader\n");
         return(1);
      }
      while(limeReaderNextRecord(limereader) != LIME_EOF )
      {
         lime_type = limeReaderType(limereader);
         if(strcmp(lime_type,"ildg-binary-data")==0)
         {
            break;
         }
         if(strcmp(lime_type,"ildg-format")==0)
         {
            lime_data_size = limeReaderBytes(limereader);
            lime_data = (char * )malloc(lime_data_size);
            limeReaderReadData((void *)lime_data,&lime_data_size, limereader);
            sscanf(qcd_getParam("<precision>",lime_data, lime_data_size),"%i",&isDouble);    
            //printf("got precision: %i\n",isDouble);
            free(lime_data);
         }
      }
      /* read 1 byte to set file-pointer to start of binary data */
      lime_data_size=1;
      limeReaderReadData(&dummy,&lime_data_size,limereader);
      offset = ftell(fid)-1;
      limeDestroyReader(limereader);      
      fclose(fid);
            
   }//end myid==0 condition
   
   MPI_Bcast(&isDouble,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&offset,sizeof(MPI_Offset),MPI_BYTE,0,MPI_COMM_WORLD);
   if(isDouble==64)
      isDouble=1;      
   else if(isDouble==32)
      isDouble=0; 
   else
   {
      fprintf(stderr,"process %i: Error in qcd_getGaugeLime! Unsupported precision\n",u->geo->myid);
      return(1);
   }   
   
   
   if(isDouble)
   {
      sizes[0]=u->geo->L[0]; 
      sizes[1]=u->geo->L[3];
      sizes[2]=u->geo->L[2];
      sizes[3]=u->geo->L[1];
      sizes[4]=4*3*3*2;
      lsizes[0]=u->geo->lL[0];
      lsizes[1]=u->geo->lL[3];
      lsizes[2]=u->geo->lL[2];
      lsizes[3]=u->geo->lL[1];
      lsizes[4]=sizes[4];
      starts[0]=u->geo->Pos[0]*lsizes[0];
      starts[1]=u->geo->Pos[3]*lsizes[1];
      starts[2]=u->geo->Pos[2]*lsizes[2];
      starts[3]=u->geo->Pos[1]*lsizes[3];                  
      starts[4]=0;

      MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock);
      MPI_Type_commit(&subblock);
      
      MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
      MPI_File_set_view(mpifid, offset, MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);
      
      //load time-slice by time-slice:
      chunksize=4*3*3*sizeof(qcd_complex_16);
      buffer = (char*) malloc(chunksize*u->geo->lV);
      if(buffer==NULL)
      {
         fprintf(stderr,"process %i: Error in qcd_getGaugeLime! Out of memory\n",u->geo->myid);
         return(1);
      }
      MPI_File_read_all(mpifid, buffer, 4*3*3*2*u->geo->lV, MPI_DOUBLE, &status);
      if(!qcd_isBigEndian())      
         qcd_swap_8((double*) buffer,2*4*3*3*u->geo->lV);
      i=0;
      for(t=0; t<u->geo->lL[0];t++)
      for(z=0; z<u->geo->lL[3];z++)
      for(y=0; y<u->geo->lL[2];y++)
      for(x=0; x<u->geo->lL[1];x++)
      for(mu=0; mu<4; mu++)
      {
         nu=(mu+1)%4;
         memcpy(&(u->D[qcd_LEXIC(t,x,y,z,u->geo->lL)][nu][0][0].re),&(buffer[i]),144);
         i+=144;
      }
      free(buffer);
      MPI_File_close(&mpifid);
      MPI_Type_free(&subblock);
      
      return 0;
   }//end isDouble condition
   else
   {
      sizes[0]=u->geo->L[0]; 
      sizes[1]=u->geo->L[3];
      sizes[2]=u->geo->L[2];
      sizes[3]=u->geo->L[1];
      sizes[4]=4*3*3*2;
      lsizes[0]=u->geo->lL[0];
      lsizes[1]=u->geo->lL[3];
      lsizes[2]=u->geo->lL[2];
      lsizes[3]=u->geo->lL[1];
      lsizes[4]=sizes[4];
      starts[0]=u->geo->Pos[0]*lsizes[0];
      starts[1]=u->geo->Pos[3]*lsizes[1];
      starts[2]=u->geo->Pos[2]*lsizes[2];
      starts[3]=u->geo->Pos[1]*lsizes[3];                  
      starts[4]=0;

      MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
      MPI_Type_commit(&subblock);
      
      MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
      MPI_File_set_view(mpifid, offset, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
      
      //load time-slice by time-slice:
      chunksize=4*3*3*2*sizeof(qcd_real_4);
      buffer = (char*) malloc(chunksize*u->geo->lV);
      if(buffer==NULL)
      {
         fprintf(stderr,"process %i: Error in qcd_getGaugeLime! Out of memory\n",u->geo->myid);
         return(1);
      }
      MPI_File_read_all(mpifid, buffer, 4*3*3*2*u->geo->lV, MPI_FLOAT, &status);
      if(!qcd_isBigEndian())      
         qcd_swap_4((float*) buffer,2*4*3*3*u->geo->lV3);
      
      for(t=0; t<u->geo->lL[0];t++)
      for(z=0; z<u->geo->lL[3];z++)
      for(y=0; y<u->geo->lL[2];y++)
      for(x=0; x<u->geo->lL[1];x++)
      for(mu=0; mu<4; mu++)
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
      {
         nu=(mu+1)%4;
         j=qcd_LEXIC(t,x,y,z,u->geo->lL);
         u->D[j][nu][c1][c2].re = *((float*) &(buffer[i])); i+=4;
         u->D[j][nu][c1][c2].im = *((float*) &(buffer[i])); i+=4;
      }      
            
      free(buffer);
      MPI_File_close(&mpifid);
      MPI_Type_free(&subblock);            

      return 0;
   }//end isDouble condition
} 
 
int qcd_getGaugeField(char *fname, int type, qcd_gaugeField *u)
{
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_getGaugeField! Gauge field must be properly initialized.\n");
      return(1);
   }
   switch(type)
   {
      case 0: // lime format
         return(qcd_getGaugeLime(fname, u));
         break;
      default:
         fprintf(stderr,"process %i: Error in qcd_getGaugeField! Unknown type.\n",u->geo->myid);
         return(1);
   }
}//end qcd_getGaugeField 




int qcd_getPropagatorCMI(char *fname, qcd_propagator *p)
{
   char vecnames[12][qcd_MAX_STRING_LENGTH];
   FILE *fid;
   qcd_int_4 i;
   qcd_uint_4 mu,nu,c1,c2;
   qcd_uint_4 x,y,z,t;
   qcd_uint_8 l;
   int sizes[5],lsizes[5],starts[5];
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   qcd_real_4 *buffer;
   
   if(!p->initialized)
   {
      fprintf(stderr,"Error in qcd_getPropagatorCMI! Propagator must be properly initialized.\n");
      return(1);
   }
   
   //open and check solution list file
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getPropagatorCMI! Cannot open %s for reading\n",p->geo->myid,fname);
      return(1);
   }
   //if(p->geo->myid==0) printf("solution list exists\n");
   i=-1;
   while(!feof(fid))
   {
      i++;
      fscanf(fid,"%s",&(vecnames[i][0]));
   }
   fclose(fid);
   if(i!=12)
   {
      fprintf(stderr,"process %i: Error in qcd_getPropagatorCMI! Expecting exactly 12 solution vectors\n",p->geo->myid);
      exit(1);
   }
   //if(p->geo->myid==0) printf("got 12 file names\n");
   for(i=0; i<12; i++)
   {
      fid=fopen(vecnames[i],"r");
      if(fid==NULL)
      {
         fprintf(stderr,"process %i: Error in qcd_getPropagatorCMI! Cannot open %s for reading\n",p->geo->myid,vecnames[i]);
         return(1);
      }
      fclose(fid);
   }
   //if(p->geo->myid==0) printf("all of them readable\n");
   
   sizes[0]=p->geo->L[1]; 
   sizes[1]=p->geo->L[2];
   sizes[2]=p->geo->L[3];
   sizes[3]=p->geo->L[0];
   sizes[4]=4*3*2;
   lsizes[0]=p->geo->lL[1];
   lsizes[1]=p->geo->lL[2];
   lsizes[2]=p->geo->lL[3];
   lsizes[3]=p->geo->lL[0];
   lsizes[4]=sizes[4];
   starts[0]=p->geo->Pos[1]*lsizes[0];
   starts[1]=p->geo->Pos[2]*lsizes[1];
   starts[2]=p->geo->Pos[3]*lsizes[2];
   starts[3]=p->geo->Pos[0]*lsizes[3];                  
   starts[4]=0;

   buffer = (qcd_real_4*) malloc(p->geo->lV*4*3*2*sizeof(qcd_real_4));
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getPropagatorCMI! Out of memory\n",p->geo->myid);
      return(1);
   }
   //if(p->geo->myid==0) printf("allocated buffer\n");
   MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
   MPI_Type_commit(&subblock);
            
   for(nu=0; nu<4; nu++)
      for(c2=0; c2<3; c2++)
      {
         MPI_File_open(MPI_COMM_WORLD, vecnames[c2+nu*3], MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
         MPI_File_set_view(mpifid, 0, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
         MPI_File_read_all(mpifid, buffer, 4*3*2*p->geo->lV, MPI_FLOAT, &status);
         if(!qcd_isBigEndian())      
            qcd_swap_4(buffer,2*4*3*p->geo->lV);
         MPI_File_close(&mpifid);         
         //reshuffle data
         l=0;
         for(x=0; x<p->geo->lL[1]; x++)
         for(y=0; y<p->geo->lL[2]; y++)
         for(z=0; z<p->geo->lL[3]; z++)
         for(t=0; t<p->geo->lL[0]; t++)
         for(mu=0; mu<4; mu++)
         for(c1=0; c1<3; c1++)
         {
            p->D[qcd_LEXIC(t,x,y,z,p->geo->lL)][mu][nu][c1][c2] = (qcd_complex_16) {buffer[l],buffer[l+1]}; l+=2;
         }
      }
   
   free(buffer);
   MPI_Type_free(&subblock);
   //if(p->geo->myid==0) printf("buffer released\n",c2+nu*3);
   return 0;
}//end qcd_getPropagatorCMI



int qcd_getPropagator(char *fname, int type, qcd_propagator *p)
{
   if(!p->initialized)
   {
      fprintf(stderr,"Error in qcd_getPropagator! Propagator must be properly initialized.\n");
      return(1);
   }
   switch(type)
   {
      case 0: // CMI format
         //if(p->geo->myid==0) printf("loading CMI type propagator\n");
         return(qcd_getPropagatorCMI(fname, p));
         break;
      default:
         fprintf(stderr,"process %i: Error in qcd_getPropagator! Unknown type.\n",p->geo->myid);
         return(1);
   }
}//end qcd_getPropagator
