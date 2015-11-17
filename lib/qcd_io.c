/* qcd_io.c
 *
 * parallel I/O routines
 * for gauge fields and 
 * propagators in several
 * formats.
 *
 * Tomasz Korzec 2008
 ******************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <mpi.h>
#include <qcd.h>
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
   return -1;
}


/*
  Function that converts between
  big endian and little endian
  8 byte words.
*/
void qcd_swap_8(double *Rd, size_t N)
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


void qcd_swap_4(float *Rd, size_t N)
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
      fgetc(pfile);
      i++;
   }
   *(len)=i;
   rewind(pfile);
   params = (char*)malloc(*(len)*sizeof(char));

   fread(params,sizeof(char),*len,pfile);
   
   fclose(pfile);
  
   return params;
}


/* by Giannis Koutsou */
char* qcd_getParam(char token[],char* params,int len)
{
   int i,token_len=strlen(token);

   for(i=0;i<len-token_len;i++)
   {
      if(memcmp(token,params+i,token_len)==0)
      {
         i+=token_len;
         *(strchr(params+i,'<'))='\0';
         break;
      }
   }
   return params+i;
}
 

/* reads the text message from the header of ILDG configuration files */
int qcd_getLimeMessage(char *fname, qcd_geometry *geo, char **message_ptr)
{
   LimeReader *limereader;
   FILE *fid;
   char *lime_type;
   n_uint64_t lime_data_size=0;
   int error_occured=0;
   char *message = NULL;
         
   if(!geo->initialized)
   {
      fprintf(stderr,"Error in qcd_getLimeMessage! Geometry must be properly initialized.\n");
      return(1);
   }   
   
   if(geo->myid == 0)
   {
      /* read lime header */
      fid=fopen(fname,"r");
      if(fid==NULL)
      {
         fprintf(stderr,"process 0: Error in qcd_getLimeMessage! Could not open %s for reading\n",fname);
         error_occured=1;
      }else
      if((limereader = limeCreateReader(fid))==NULL)
      {
         fprintf(stderr,"process 0: Error in qcd_getLimeMessage! Could not create limeReader\n");
         error_occured=1;
      }
      
      if(!error_occured)
      {
         while(limeReaderNextRecord(limereader) != LIME_EOF)
         {
            lime_type = limeReaderType(limereader);
            if(strcmp(lime_type,"ildg-binary-data")==0)
            {
               break;
            }
            if(strcmp(lime_type,"xlf-info")==0)
            {
               lime_data_size = limeReaderBytes(limereader);
               message = (char * )malloc(lime_data_size+1);
               limeReaderReadData((void *)message,&lime_data_size, limereader);
               message[lime_data_size]='\0';
               //printf("lime-message:%s",message);
            }
         }
         fclose(fid);
      }
   }//end myid==0 condition
   MPI_Bcast(&error_occured,1,MPI_INT,0,MPI_COMM_WORLD);
   if(error_occured) return(1);
   
   MPI_Bcast(&lime_data_size,sizeof(n_uint64_t),MPI_BYTE,0,MPI_COMM_WORLD);
   if(lime_data_size>0)
   {
      if(geo->myid != 0)
         message = (char *) malloc(lime_data_size+1);
      MPI_Bcast(&(message[0]),lime_data_size+1,MPI_BYTE,0,MPI_COMM_WORLD);
      //printf("process %i: lime-message:%s\n",geo->myid,message);
   }else
   {
      message = (char *) malloc(1);
      message[0]='\0';
   }   
   *message_ptr = message;
   return(0);
}//end qcd_getLimeGmessage


 
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
   qcd_uint_8 i=0,j=0;
   qcd_uint_2 chunksize,mu,nu,c1,c2;
   char *buffer;
   qcd_uint_4 x,y,z,t;
   int  isDouble;
   int error_occured=0;
      
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
         error_occured=1;
      }
      if ((limereader = limeCreateReader(fid))==NULL)
      {
         fprintf(stderr,"process 0: Error in qcd_getGaugeLime! Could not create limeReader\n");
         error_occured=1;
      }
      if(!error_occured)
      {
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
      }     
   }//end myid==0 condition
   
   MPI_Bcast(&error_occured,1,MPI_INT,0,MPI_COMM_WORLD);
   if(error_occured) return(1);
   
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
	qcd_swap_8((double*) buffer,(size_t)(2*4*3*3)*(size_t)u->geo->lV);
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
	qcd_swap_4((float*) buffer,(size_t)(2*4*3*3)*(size_t)u->geo->lV);
      
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
}//end qcd_getGaugeLime 
 
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


/*writes a double precision ildg formatted gauge configuration file*/
int qcd_writeGaugeLime(char *fname, qcd_gaugeField *u, char *message)
{
   FILE *fid;   
   int error_in_header=0;
   LimeWriter *limewriter;
   LimeRecordHeader *limeheader;
   int ME_flag=0, MB_flag=0, limeStatus;
   qcd_uint_8 message_length;
   MPI_Offset offset;
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   int sizes[5], lsizes[5], starts[5];
   qcd_uint_8 i;
   qcd_uint_2 chunksize,mu,nu;
   char *buffer;
   qcd_uint_4 x,y,z,t;
   char ildgHeader[2048];

   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_getGaugeLime! Gauge field must be properly initialized.\n");
      return(1);
   }   


   //process 0 creates file and writes header   
   if(u->geo->myid == 0)
   {
      fid=fopen(fname,"w");
      if(fid==NULL)
      {
         fprintf(stderr,"process 0: Error in qcd_writeGaugeLime! Could not open %s for writing\n",fname);
         error_in_header=1;
      }else
      {

         limewriter = limeCreateWriter(fid); 
         if(limewriter == (LimeWriter*)NULL) {
            fprintf(stderr, "Error in qcd_writeGaugeLime. LIME error in file %s for writing!\n", fname);
            error_in_header=1;
         }else
         {            
            sprintf(ildgHeader, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version> 1.0 </version>\n  <field> su3gauge </field>\n  <precision> 64 </precision>\n  <lx> %d </lx>\n  <ly> %d </ly>\n  <lz> %d </lz>\n  <lt> %d </lt>\n</ildgFormat>",u->geo->L[1],u->geo->L[2], u->geo->L[3], u->geo->L[0]);
            message_length=(qcd_uint_8) strlen(ildgHeader);
            MB_flag=1; ME_flag=0;
            limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-format", message_length);
            if(limeheader == (LimeRecordHeader*)NULL) 
            {
               fprintf(stderr, "Error in qcd_writeGaugeLime. LIME create header ildg-format error.\n");
               error_in_header=1;
            }
            limeStatus = limeWriteRecordHeader(limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in qcd_writeGaugeLime. LIME write header ildg-format error %d\n", limeStatus);
              error_in_header=1;
            }
            limeDestroyHeader(limeheader);
            limeStatus = limeWriteRecordData(ildgHeader, &message_length, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in qcd_writeGaugeLime. LIME write header ildg-format error %d\n", limeStatus);
              error_in_header=1;
            }
            
            message_length=strlen(message);
            ME_flag=0; MB_flag=0;
            limeheader = limeCreateHeader(MB_flag, ME_flag, "xlf-info", message_length);
            limeStatus = limeWriteRecordHeader(limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
               fprintf(stderr, "Error in qcd_writeGaugeLime. LIME write header error %d\n", limeStatus);
               error_in_header=1;
            }
            limeDestroyHeader( limeheader );
            limeWriteRecordData(message, &message_length, limewriter);

            message_length = u->geo->V*4*3*3*sizeof(qcd_complex_16);
            MB_flag=0; ME_flag=1;
            limeheader = limeCreateHeader(MB_flag, ME_flag, "ildg-binary-data", message_length);
            limeStatus = limeWriteRecordHeader( limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
               fprintf(stderr, "Error in qcd_writeGaugeLime. LIME write header error %d\n", limeStatus);
               error_in_header=1;
            }
            limeDestroyHeader( limeheader );
         }
         //make one fake record-write to set offset
         message_length=1;
         limeWriteRecordData(message, &message_length, limewriter);
         offset = ftell(fid)-1;
         fclose(fid);
      }//fid non null condition
   }//myid==0 condition
   
   MPI_Bcast(&error_in_header,1,MPI_INT,0,MPI_COMM_WORLD);
   if(error_in_header)
      return(1);
   MPI_Bcast(&offset,sizeof(MPI_Offset),MPI_BYTE,0,MPI_COMM_WORLD);
   
   //header written. Now use MPI-I/O to write binary data
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

   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, offset, MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);

   chunksize=4*3*3*sizeof(qcd_complex_16);
   buffer = (char*) malloc(chunksize*u->geo->lV);
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_writeGaugeLime! Out of memory\n",u->geo->myid);
      return(1);
   }
   i=0;
   for(t=0; t<u->geo->lL[0];t++)
   for(z=0; z<u->geo->lL[3];z++)
   for(y=0; y<u->geo->lL[2];y++)
   for(x=0; x<u->geo->lL[1];x++)
   for(mu=0; mu<4; mu++)
   {
      nu=(mu+1)%4;
      memcpy(&(buffer[i]),&(u->D[qcd_LEXIC(t,x,y,z,u->geo->lL)][nu][0][0].re),144);
      i+=144;
   }
   if(!qcd_isBigEndian())
     qcd_swap_8((double*) buffer,(size_t)(2*4*3*3)*(size_t)u->geo->lV);
   
   MPI_File_write_all(mpifid, buffer, 4*3*3*2*u->geo->lV, MPI_DOUBLE, &status);

   free(buffer);
   MPI_File_close(&mpifid);
   MPI_Type_free(&subblock);

   return 0;
}//end qcd_writeGaugeLime


int qcd_writeGaugeField(char *fname, int type, qcd_gaugeField *u, ...)
{
   va_list arg_list;
   char *message = NULL;
   
   va_start(arg_list, u);
   if(type==0)
      message = va_arg(arg_list,char *);
   va_end(arg_list);   
   
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_writeGaugeField! Gauge field must be properly initialized.\n");
      return(1);
   }
   switch(type)
   {
      case 0: // lime format
         return(qcd_writeGaugeLime(fname, u,message));
         break;
      default:
         fprintf(stderr,"process %i: Error in qcd_writeGaugeField! Unknown type.\n",u->geo->myid);
         return(1);
   }
}//end qcd_writeGaugeField 



int qcd_getPropagatorCMI(char *fname, qcd_propagator *p)
{
   char vecnames[12][qcd_MAX_STRING_LENGTH];
   FILE *fid;
   qcd_uint_4 mu,nu,c1,c2,i;
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
	   qcd_swap_4(buffer,(size_t)(2*4*3)*(size_t)p->geo->lV);
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

int qcd_getPropagatorLime(char *fname, qcd_propagator *p)
{
   char vecnames[12][qcd_MAX_STRING_LENGTH];
   FILE *fid;
   qcd_int_4 i;
   qcd_uint_4 mu,nu,c1,c2;
   qcd_uint_4 x,y,z,t;
   qcd_vector *v = malloc (sizeof (qcd_vector));
   
   if(!p->initialized)
   {
     fprintf(stderr,"Error in %s! Propagator must be properly initialized.\n", __func__);
      return(1);
   }
   
   //open and check solution list file
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
     fprintf(stderr,"process %i: Error in %s! Cannot open %s for reading\n",p->geo->myid,__func__, fname);
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
     fprintf(stderr,"process %i: Error in %s! Expecting exactly 12 solution vectors\n",p->geo->myid, __func__);
      exit(1);
   }
   //if(p->geo->myid==0) printf("got 12 file names\n");
   for(nu=0; nu<4; nu++)
   for(c2=0; c2<3; c2++)
     {
       i=c2 + nu*3;
      fid=fopen(vecnames[i],"r");
      if(fid==NULL)
      {
         fprintf(stderr,"process %i: Error in %s! Cannot open %s for reading\n",p->geo->myid, __func__,vecnames[i]);
         return(1);
      }
      fclose(fid);

      qcd_initVector (v, p->geo);

      qcd_getVectorLime (vecnames[i], v);
      for(x=0; x<p->geo->lL[1]; x++)
	for(y=0; y<p->geo->lL[2]; y++)
	  for(z=0; z<p->geo->lL[3]; z++)
	    for(t=0; t<p->geo->lL[0]; t++)
	      for(mu=0; mu<4; mu++)
		for(c1=0; c1<3; c1++)
		  {
		    p->D[qcd_LEXIC(t,x,y,z,p->geo->lL)][mu][nu][c1][c2] = 
		      (qcd_complex_16) {v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1].re,
					v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1].im};
		  }
      qcd_destroyVector (v);
     }
      
   free (v);
   return 0;
}//end qcd_getPropagatorLime

int qcd_getPropagatorHMC(char *fname, qcd_propagator *p)
{
   FILE *fid;
   double beta,kappa,mu_tm;
   int fL, fT;
   long firstlinelength;
   qcd_uint_4 mu,nu,c1,c2;
   qcd_uint_4 x,y,z,t;
   qcd_uint_8 l;
   int sizes[6],lsizes[6],starts[6];
   MPI_Datatype subblock6d;  //MPI-type, 6d subarray
   MPI_File mpifid;
   MPI_Status status;
   qcd_real_8 *buffer;
   
   if(!p->initialized)
   {
      fprintf(stderr,"Error in qcd_getPropagatorHMC! Propagator must be properly initialized.\n");
      return(1);
   }
   
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getPropagatorHMC! Cannot open %s for reading\n",p->geo->myid,fname);
      return(1);
   }
   fscanf(fid,"%lf %lf %lf %d %d\n", &beta, &kappa, &mu_tm, &fL, &fT);
   if((fL!=p->geo->L[1])||(fT!=p->geo->L[0]))
   {
      fprintf(stderr, "process %i: Error in qcd_getPropagatorHMC! Dimensions in file %s don't match geometry!\n",p->geo->myid, fname);
      return(1);
   }
   firstlinelength=ftell(fid);
   //if(p->geo->myid==0) printf("firstlinelength = %i\n",firstlinelength);
   fclose(fid);
   
   sizes[0]=4*3;
   sizes[1]=p->geo->L[1]; 
   sizes[2]=p->geo->L[2];
   sizes[3]=p->geo->L[3];
   sizes[4]=p->geo->L[0];
   sizes[5]=4*3*2;
   lsizes[0]=4*3;
   lsizes[1]=p->geo->lL[1];
   lsizes[2]=p->geo->lL[2];
   lsizes[3]=p->geo->lL[3];
   lsizes[4]=p->geo->lL[0];
   lsizes[5]=4*3*2;
   starts[0]=0;
   starts[1]=p->geo->Pos[1]*lsizes[1];
   starts[2]=p->geo->Pos[2]*lsizes[2];
   starts[3]=p->geo->Pos[3]*lsizes[3];
   starts[4]=p->geo->Pos[0]*lsizes[4];                  
   starts[5]=0;

   buffer = (qcd_real_8*) malloc(4*3*p->geo->lV*4*3*2*sizeof(qcd_real_8));
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getPropagatorHMC! Out of memory\n",p->geo->myid);
      return(1);
   }
   //if(p->geo->myid==0) printf("allocated buffer\n");
   MPI_Type_create_subarray(6,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock6d);
   MPI_Type_commit(&subblock6d);

   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, (MPI_Offset) firstlinelength, MPI_DOUBLE, subblock6d, "native", MPI_INFO_NULL);

   MPI_File_read_all(mpifid, buffer, 4*3*4*3*2*p->geo->lV, MPI_DOUBLE, &status);
   if(!qcd_isBigEndian())      
     qcd_swap_8(buffer,(size_t)(4*3*4*3*2)*(size_t)p->geo->lV);
   l=0;
   //reshuffle data
   for(nu=0; nu<4; nu++)
   for(c2=0; c2<3; c2++)
   for(x=0; x<p->geo->lL[1]; x++)
   for(y=0; y<p->geo->lL[2]; y++)
   for(z=0; z<p->geo->lL[3]; z++)
   for(t=0; t<p->geo->lL[0]; t++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   {
      p->D[qcd_LEXIC(t,x,y,z,p->geo->lL)][mu][nu][c1][c2] = (qcd_complex_16) {buffer[l],buffer[l+1]}; l+=2;
   }
   
   
   MPI_File_close(&mpifid);            
   
   free(buffer);
   MPI_Type_free(&subblock6d);
   //if(p->geo->myid==0) printf("buffer released\n",c2+nu*3);
   return 0;
}//end qcd_getPropagatorHMC



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
      case 1: // old HMC format
         //if(p->geo->myid==0) printf("loading HMC type propagator\n");
         return(qcd_getPropagatorHMC(fname, p));
         break;   
      case qcd_PROP_LIME:
         return(qcd_getPropagatorLime(fname, p));
         break;   
      default:
         fprintf(stderr,"process %i: Error in qcd_getPropagator! Unknown type.\n",p->geo->myid);
         return(1);
   }
}//end qcd_getPropagator



int qcd_getVectorCMI(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v)
{
   FILE *fid;
   qcd_uint_2 mu,c1;
   qcd_uint_4 x,y,z,t;
   qcd_uint_8 l;
   int sizes[5],lsizes[5],starts[5];
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   qcd_real_4 *buffer;
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_getVectorCMI! Vector must be properly initialized.\n");
      return(1);
   }
   
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getVectorCMI! Cannot open %s for reading\n",v->geo->myid,fname);
      return(1);
   }
   fclose(fid);
   
   sizes[0]=v->geo->L[1]; 
   sizes[1]=v->geo->L[2];
   sizes[2]=v->geo->L[3];
   sizes[3]=v->geo->L[0];
   sizes[4]=4*3*2;
   lsizes[0]=v->geo->lL[1];
   lsizes[1]=v->geo->lL[2];
   lsizes[2]=v->geo->lL[3];
   lsizes[3]=v->geo->lL[0];
   lsizes[4]=sizes[4];
   starts[0]=v->geo->Pos[1]*lsizes[0];
   starts[1]=v->geo->Pos[2]*lsizes[1];
   starts[2]=v->geo->Pos[3]*lsizes[2];
   starts[3]=v->geo->Pos[0]*lsizes[3];                  
   starts[4]=0;

   buffer = (qcd_real_4*) malloc(v->geo->lV*4*3*2*sizeof(qcd_real_4));
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getVectorCMI! Out of memory\n",v->geo->myid);
      return(1);
   }
   MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
   MPI_Type_commit(&subblock);
            
   
   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, 0, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
   MPI_File_read_all(mpifid, buffer, 4*3*2*v->geo->lV, MPI_FLOAT, &status);
   if(!qcd_isBigEndian())      
     qcd_swap_4(buffer,(size_t)(2*4*3)*(size_t)v->geo->lV);
   MPI_File_close(&mpifid);
   //reshuffle data
   l=0;
   for(x=0; x<v->geo->lL[1]; x++)
   for(y=0; y<v->geo->lL[2]; y++)
   for(z=0; z<v->geo->lL[3]; z++)
   for(t=0; t<v->geo->lL[0]; t++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   {
      v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1] = (qcd_complex_16) {buffer[l],buffer[l+1]}; l+=2;
   }

   
   free(buffer);
   MPI_Type_free(&subblock);
   return 0;
}//end qcd_getVectorCMI

int qcd_getVectorLime(char *fname, qcd_vector *v)
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
   qcd_uint_2 chunksize,mu,c1;
   char *buffer;
   qcd_uint_4 x,y,z,t;
   int  isDouble;
   int error_occured=0;
   int next_rec_is_prop = 0;

   if(!v->initialized)
   {
     fprintf(stderr,"Error in %s! Gauge field must be properly initialized.\n", __func__);
      return(1);
   }   
   
   if(v->geo->myid == 0)
   {
      /* read lime header */
      fid=fopen(fname,"r");
      if(fid==NULL)
      {
	fprintf(stderr,"process 0: Error in %s Could not open %s for reading\n",__func__, fname);
         error_occured=1;
      }
      if ((limereader = limeCreateReader(fid))==NULL)
      {
	fprintf(stderr,"process 0: Error in %s! Could not create limeReader\n", __func__);
         error_occured=1;
      }
      if(!error_occured)
      {
	while(limeReaderNextRecord(limereader) != LIME_EOF )
         {
	   lime_type = limeReaderType(limereader);
            if(strcmp(lime_type,"propagator-type")==0)
            {
	      lime_data_size = limeReaderBytes(limereader);
	      lime_data = (char * )malloc(lime_data_size);
	      limeReaderReadData((void *)lime_data,&lime_data_size, limereader);

	      if (strncmp ("DiracFermion_Source_Sink_Pairs", lime_data, strlen ("DiracFermion_Source_Sink_Pairs"))!=0 &&
		  strncmp ("DiracFermion_Sink", lime_data, strlen ("DiracFermion_Sink"))!=0 )
		{
		  fprintf (stderr, " process 0: Error in %s! Got %s for \"propagator-type\", expecting %s or %s\n", __func__, lime_data, 
			   "DiracFermion_Source_Sink_Pairs", 
			   "DiracFermion_Sink");
		  error_occured = 1;
		  break;
		}
	      free(lime_data);
            }
//lime_type="scidac-binary-data";
            if((strcmp(lime_type,"etmc-propagator-format")==0) || (strcmp(lime_type,"etmc-source-format")==0))
            {
	      lime_data_size = limeReaderBytes(limereader);
	      lime_data = (char * )malloc(lime_data_size);
	      limeReaderReadData((void *)lime_data,&lime_data_size, limereader);
	      sscanf(qcd_getParam("<precision>",lime_data, lime_data_size),"%i",&isDouble);    
	      //printf("got precision: %i\n",isDouble);
	      free(lime_data);
	      
	      next_rec_is_prop = 1;
            }
            if(strcmp(lime_type,"scidac-binary-data")==0 && next_rec_is_prop)
            {	      
	      break;
            }
         }
         /* read 1 byte to set file-pointer to start of binary data */
         lime_data_size=1;
         limeReaderReadData(&dummy,&lime_data_size,limereader);
         offset = ftell(fid)-1;
         limeDestroyReader(limereader);      
         fclose(fid);
      }     
   }//end myid==0 condition
   
   MPI_Bcast(&error_occured,1,MPI_INT,0,MPI_COMM_WORLD);
   if(error_occured) return(1);
   
   MPI_Bcast(&isDouble,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Bcast(&offset,sizeof(MPI_Offset),MPI_BYTE,0,MPI_COMM_WORLD);
   if(isDouble==64)
      isDouble=1;      
   else if(isDouble==32)
      isDouble=0; 
   else
   {
     fprintf(stderr,"process %i: Error in %s! Unsupported precision\n",v->geo->myid, __func__);
   }  
   
   if(isDouble)
   {
      sizes[0]=v->geo->L[0]; 
      sizes[1]=v->geo->L[3];
      sizes[2]=v->geo->L[2];
      sizes[3]=v->geo->L[1];
      sizes[4]=4*3*2;
      lsizes[0]=v->geo->lL[0];
      lsizes[1]=v->geo->lL[3];
      lsizes[2]=v->geo->lL[2];
      lsizes[3]=v->geo->lL[1];
      lsizes[4]=sizes[4];
      starts[0]=v->geo->Pos[0]*lsizes[0];
      starts[1]=v->geo->Pos[3]*lsizes[1];
      starts[2]=v->geo->Pos[2]*lsizes[2];
      starts[3]=v->geo->Pos[1]*lsizes[3];                  
      starts[4]=0;

      MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock);
      MPI_Type_commit(&subblock);
      
      MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
      MPI_File_set_view(mpifid, offset, MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);
      
      //load time-slice by time-slice:
      chunksize=4*3*2*sizeof(qcd_real_8);
      buffer = (char*) malloc(chunksize*v->geo->lV);
      if(buffer==NULL)
      {
	fprintf(stderr,"process %i: Error in %s! Out of memory\n",v->geo->myid, __func__);
         return(1);
      }
      MPI_File_read_all(mpifid, buffer, 4*3*2*v->geo->lV, MPI_DOUBLE, &status);
      if(!qcd_isBigEndian())      
	qcd_swap_8((double*) buffer,(size_t)(2*4*3)*(size_t)v->geo->lV);
      i=0;
      for(t=0; t<v->geo->lL[0];t++)
      for(z=0; z<v->geo->lL[3];z++)
      for(y=0; y<v->geo->lL[2];y++)
      for(x=0; x<v->geo->lL[1];x++)
	for(mu=0; mu<4; mu++)
	for(c1=0; c1<3; c1++)
	{
	  v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1] = (qcd_complex_16){
	    ((qcd_real_8 *)buffer)[i], 
	    ((qcd_real_8 *)buffer)[i+1]};
	  i+=2;
      }
      free(buffer);
      MPI_File_close(&mpifid);
      MPI_Type_free(&subblock);
      
      return 0;
   }//end isDouble condition
   else
   {
      sizes[0]=v->geo->L[0]; 
      sizes[1]=v->geo->L[3];
      sizes[2]=v->geo->L[2];
      sizes[3]=v->geo->L[1];
      sizes[4]=4*3*2;
      lsizes[0]=v->geo->lL[0];
      lsizes[1]=v->geo->lL[3];
      lsizes[2]=v->geo->lL[2];
      lsizes[3]=v->geo->lL[1];
      lsizes[4]=sizes[4];
      starts[0]=v->geo->Pos[0]*lsizes[0];
      starts[1]=v->geo->Pos[3]*lsizes[1];
      starts[2]=v->geo->Pos[2]*lsizes[2];
      starts[3]=v->geo->Pos[1]*lsizes[3];                  
      starts[4]=0;

      MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
      MPI_Type_commit(&subblock);
      
      MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
      MPI_File_set_view(mpifid, offset, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
      
      //load time-slice by time-slice:
      chunksize=4*3*2*sizeof(qcd_real_4);
      buffer = (char*) malloc(chunksize*v->geo->lV);
      if(buffer==NULL)
      {
	fprintf(stderr,"process %i: Error in %s! Out of memory\n",v->geo->myid, __func__);
         return(1);
      }
      MPI_File_read_all(mpifid, buffer, 4*3*2*v->geo->lV, MPI_FLOAT, &status);

      if(!qcd_isBigEndian())
	qcd_swap_4((float*) buffer,(size_t)(2*4*3)*(size_t)v->geo->lV);

      i=0;
      for(t=0; t<v->geo->lL[0];t++)
      for(z=0; z<v->geo->lL[3];z++)
      for(y=0; y<v->geo->lL[2];y++)
      for(x=0; x<v->geo->lL[1];x++)
      for(mu=0; mu<4; mu++)
      for(c1=0; c1<3; c1++)
      {
	j=qcd_LEXIC(t,x,y,z,v->geo->lL);
	v->D[j][mu][c1].re = *((float*) (buffer+i)); i+=4;
	v->D[j][mu][c1].im = *((float*) (buffer+i)); i+=4;
      }      
            
      free(buffer);
      MPI_File_close(&mpifid);
      MPI_Type_free(&subblock);            

      return 0;
   }//end isDouble condition
}//end qcd_getVectorLime 


int qcd_getVectorHMC(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v)
{
   FILE *fid;
   double beta,kappa,mu_tm;
   int fL, fT;
   long firstlinelength;
   qcd_uint_2 mu,c1;
   qcd_uint_4 x,y,z,t;
   qcd_uint_8 l;
   int sizes[6],lsizes[6],starts[6];
   MPI_Datatype subblock6d;  //MPI-type, 6d subarray
   MPI_File mpifid;
   MPI_Status status;
   qcd_real_8 *buffer;
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_getVectorHMC! Vector must be properly initialized.\n");
      return(1);
   }
   
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getVectorHMC! Cannot open %s for reading\n",v->geo->myid,fname);
      return(1);
   }
   fscanf(fid,"%lf %lf %lf %d %d\n", &beta, &kappa, &mu_tm, &fL, &fT);
   if((fL!=v->geo->L[1])||(fT!=v->geo->L[0]))
   {
      fprintf(stderr, "process %i: Error in qcd_getVectorHMC! Dimensions in file %s don't match geometry!\n",v->geo->myid, fname);
      return(1);
   }
   firstlinelength=ftell(fid);
   fclose(fid);
   
   sizes[0]=4*3;
   sizes[1]=v->geo->L[1]; 
   sizes[2]=v->geo->L[2];
   sizes[3]=v->geo->L[3];
   sizes[4]=v->geo->L[0];
   sizes[5]=4*3*2;
   lsizes[0]=4*3;
   lsizes[1]=v->geo->lL[1];
   lsizes[2]=v->geo->lL[2];
   lsizes[3]=v->geo->lL[3];
   lsizes[4]=v->geo->lL[0];
   lsizes[5]=4*3*2;
   starts[0]=0;
   starts[1]=v->geo->Pos[1]*lsizes[1];
   starts[2]=v->geo->Pos[2]*lsizes[2];
   starts[3]=v->geo->Pos[3]*lsizes[3];
   starts[4]=v->geo->Pos[0]*lsizes[4];                  
   starts[5]=0;

   buffer = (qcd_real_8*) malloc(v->geo->lV*4*3*2*sizeof(qcd_real_8));
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getVectorHMC! Out of memory\n",v->geo->myid);
      return(1);
   }
   MPI_Type_create_subarray(6,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock6d);
   MPI_Type_commit(&subblock6d);

   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, (MPI_Offset) firstlinelength, MPI_DOUBLE, subblock6d, "native", MPI_INFO_NULL);


   MPI_File_read_at_all(mpifid, (MPI_Offset) (nu*3+c2)*4*3*2*v->geo->lV, buffer, 4*3*2*v->geo->lV, MPI_DOUBLE, &status);
   if(!qcd_isBigEndian())      
     qcd_swap_8(buffer,(size_t)(4*3*2)*(size_t)v->geo->lV);
   l=0;
   //reshuffle data
   for(x=0; x<v->geo->lL[1]; x++)
   for(y=0; y<v->geo->lL[2]; y++)
   for(z=0; z<v->geo->lL[3]; z++)
   for(t=0; t<v->geo->lL[0]; t++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   {
      v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1] = (qcd_complex_16) {buffer[l],buffer[l+1]}; l+=2;
   }
   
   
   MPI_File_close(&mpifid);            
   
   free(buffer);
   MPI_Type_free(&subblock6d);
   return 0;
}//end qcd_getVectorHMC


int qcd_getVectorHMCV(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v)
{
   FILE *fid;
   double beta,kappa,mu_tm;
   int fL, fT;
   long firstlinelength;
   qcd_uint_2 mu,c1;
   qcd_uint_4 x,y,z,t;
   qcd_uint_8 l;
   int sizes[6],lsizes[6],starts[6];
   MPI_Datatype subblock6d;  //MPI-type, 6d subarray
   MPI_File mpifid;
   MPI_Status status;
   qcd_real_8 *buffer;
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_getVectorHMCV! Vector must be properly initialized.\n");
      return(1);
   }
   
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getVectorHMCV! Cannot open %s for reading\n",v->geo->myid,fname);
      return(1);
   }
   fscanf(fid,"%lf %lf %lf %d %d\n", &beta, &kappa, &mu_tm, &fL, &fT);
   if((fL!=v->geo->L[1])||(fT!=v->geo->L[0]))
   {
      fprintf(stderr, "process %i: Error in qcd_getVectorHMCV! Dimensions in file %s don't match geometry!\n",v->geo->myid, fname);
      return(1);
   }
   firstlinelength=ftell(fid);
   fclose(fid);
   
   sizes[0]=4*3;
   sizes[1]=v->geo->L[1]; 
   sizes[2]=v->geo->L[2];
   sizes[3]=v->geo->L[3];
   sizes[4]=v->geo->L[0];
   sizes[5]=4*3*2;
   lsizes[0]=4*3;
   lsizes[1]=v->geo->lL[1];
   lsizes[2]=v->geo->lL[2];
   lsizes[3]=v->geo->lL[3];
   lsizes[4]=v->geo->lL[0];
   lsizes[5]=4*3*2;
   starts[0]=0;
   starts[1]=v->geo->Pos[1]*lsizes[1];
   starts[2]=v->geo->Pos[2]*lsizes[2];
   starts[3]=v->geo->Pos[3]*lsizes[3];
   starts[4]=v->geo->Pos[0]*lsizes[4];                  
   starts[5]=0;

   buffer = (qcd_real_8*) malloc(v->geo->lV*4*3*2*sizeof(qcd_real_8));
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_getVectorHMCV! Out of memory\n",v->geo->myid);
      return(1);
   }
   MPI_Type_create_subarray(6,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock6d);
   MPI_Type_commit(&subblock6d);

   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, (MPI_Offset) firstlinelength, MPI_DOUBLE, subblock6d, "native", MPI_INFO_NULL);


   MPI_File_read_at_all(mpifid, (MPI_Offset) 0, buffer, 4*3*2*v->geo->lV, MPI_DOUBLE, &status);
   if(!qcd_isBigEndian())      
     qcd_swap_8(buffer,(size_t)(4*3*2)*(size_t)v->geo->lV);
   l=0;
   //reshuffle data
   for(x=0; x<v->geo->lL[1]; x++)
   for(y=0; y<v->geo->lL[2]; y++)
   for(z=0; z<v->geo->lL[3]; z++)
   for(t=0; t<v->geo->lL[0]; t++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   {
      v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1] = (qcd_complex_16) {buffer[l],buffer[l+1]}; l+=2;
   }
   
   
   MPI_File_close(&mpifid);            
   
   free(buffer);
   MPI_Type_free(&subblock6d);
   return 0;
}//end qcd_getVectorHMCV


int qcd_getVector(char *fname, int type, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v)
{
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_getVector! Vector must be properly initialized.\n");
      return(1);
   }
   switch(type)
   {
      case 0: // CMI format
         return(qcd_getVectorCMI(fname, nu, c2, v));
         break;
      case 1: // old HMC format
         return(qcd_getVectorHMC(fname, nu, c2, v));
         break;  
      case 2: // old HMC format with 1 file/vector 
         return(qcd_getVectorHMCV(fname, nu, c2, v));
         break;
   case qcd_PROP_LIME:
     return(qcd_getVectorLime(fname, v));
         break;
      default:
         fprintf(stderr,"process %i: Error in qcd_getVector! Unknown type.\n",v->geo->myid);
         return(1);
   }
}//end qcd_getVector




int qcd_writePropagatorCMI(char *fname, qcd_propagator *p)
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
      fprintf(stderr,"Error in qcd_writePropagatorCMI! Propagator must be properly initialized.\n");
      return(1);
   }
   
   //open and check solution list file
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_writePropagatorCMI! Cannot open %s for reading\n",p->geo->myid,fname);
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
      fprintf(stderr,"process %i: Error in qcd_writePropagatorCMI! Expecting exactly 12 vectors\n",p->geo->myid);
      exit(1);
   }
   //if(p->geo->myid==0) printf("got 12 file names\n");
   for(i=0; i<12; i++)
   {
      fid=fopen(vecnames[i],"w");
      if(fid==NULL)
      {
         fprintf(stderr,"process %i: Error in qcd_writePropagatorCMI! Cannot open %s for writing\n",p->geo->myid,vecnames[i]);
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
      fprintf(stderr,"process %i: Error in qcd_writePropagatorCMI! Out of memory\n",p->geo->myid);
      return(1);
   }
   //if(p->geo->myid==0) printf("allocated buffer\n");
   MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
   MPI_Type_commit(&subblock);
            
   for(nu=0; nu<4; nu++)
      for(c2=0; c2<3; c2++)
      {
         MPI_File_open(MPI_COMM_WORLD, vecnames[c2+nu*3], MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifid);
         MPI_File_set_view(mpifid, 0, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
         //reshuffle data
         l=0;
         for(x=0; x<p->geo->lL[1]; x++)
         for(y=0; y<p->geo->lL[2]; y++)
         for(z=0; z<p->geo->lL[3]; z++)
         for(t=0; t<p->geo->lL[0]; t++)
         for(mu=0; mu<4; mu++)
         for(c1=0; c1<3; c1++)
         {
            buffer[l++] = p->D[qcd_LEXIC(t,x,y,z,p->geo->lL)][mu][nu][c1][c2].re;
            buffer[l++] = p->D[qcd_LEXIC(t,x,y,z,p->geo->lL)][mu][nu][c1][c2].im;
         }
         if(!qcd_isBigEndian())      
	   qcd_swap_4(buffer,(size_t)(2*4*3)*(size_t)p->geo->lV);
         
         MPI_File_write_all(mpifid, buffer, 4*3*2*p->geo->lV, MPI_FLOAT, &status);
         MPI_File_close(&mpifid);         
      }
   
   free(buffer);
   MPI_Type_free(&subblock);
   //if(p->geo->myid==0) printf("buffer released\n",c2+nu*3);
   return 0;
}//end qcd_writePropagatorCMI


int qcd_writePropagatorHMCV(char *fname, qcd_propagator *p)
{
   char vecnames[12][qcd_MAX_STRING_LENGTH];
   char headerline[qcd_MAX_STRING_LENGTH];
   FILE *fid;
   long firstlinelength;
   qcd_int_4 i;
   qcd_uint_4 mu,nu,c1,c2;
   qcd_uint_4 x,y,z,t;
   qcd_uint_8 l;
   int sizes[5],lsizes[5],starts[5];
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   qcd_real_8 *buffer;
   
   if(!p->initialized)
   {
      fprintf(stderr,"Error in qcd_writePropagatorHMCV! Propagator must be properly initialized.\n");
      return(1);
   }
   
   //open and check solution list file
   fid=fopen(fname,"r");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_writePropagatorHMCV! Cannot open %s for reading\n",p->geo->myid,fname);
      return(1);
   }
//   if(p->geo->myid==0) printf("solution list exists\n");
   i=-1;
   fgets(headerline, qcd_MAX_STRING_LENGTH, fid);
//   if(p->geo->myid==0) printf("got headerline: %s",headerline);
   while(!feof(fid))
   {
      i++;
      fscanf(fid,"%s\n",&(vecnames[i][0]));
//      if(p->geo->myid==0) printf("i=%i, vecnames[i] = %s\n",i,vecnames[i]);
   }
//   if(p->geo->myid==0) printf("read %i vector file names\n",i);
   fclose(fid);
   if(i!=11)
   {
      fprintf(stderr,"process %i: Error in qcd_writePropagatorHMCV! Expecting exactly 12 vectors\n",p->geo->myid);
      exit(1);
   }
   if(p->geo->myid==0) printf("got 12 file names\n");
/* 
  for(i=0; i<12; i++)
   {
      fid=fopen(vecnames[i],"w");
      if(fid==NULL)
      {
         fprintf(stderr,"process %i: Error in qcd_writePropagatorHMCV! Cannot open %s for writing\n",p->geo->myid,vecnames[i]);
         return(1);
      }
      if(p->geo->myid==0)
      {
         fprintf(fid,"%s", headerline);
         firstlinelength=ftell(fid);
      }
      fclose(fid);      
   }
*/
   firstlinelength = strlen(headerline);
//   MPI_Bcast(&firstlinelength, 1, MPI_LONG, 0, MPI_COMM_WORLD);
//   printf("process %i: firstlinelength=%i\n",p->geo->myid,firstlinelength);  
 
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

   buffer = (qcd_real_8*) malloc(p->geo->lV*4*3*2*sizeof(qcd_real_8));
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_writePropagatorHMCV! Out of memory\n",p->geo->myid);
      return(1);
   }
//   if(p->geo->myid==0) printf("allocated buffer\n");
   MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock);
   MPI_Type_commit(&subblock);
            
   for(nu=0; nu<4; nu++)
      for(c2=0; c2<3; c2++)
      {
         MPI_File_open(MPI_COMM_WORLD, vecnames[c2+nu*3], MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &mpifid);
         if(mpifid==NULL)
         {
            fprintf(stderr,"process %i: Error in qcd_writePropagatorHMCV! Could not open output file\n",p->geo->myid);
            return(1);
         }
         if(p->geo->myid==0)
         {
            MPI_File_write(mpifid, headerline, firstlinelength, MPI_CHAR, &status);
         }
         MPI_Barrier(MPI_COMM_WORLD);
         MPI_File_sync(mpifid);
         MPI_File_set_view(mpifid, (MPI_Offset) firstlinelength, MPI_DOUBLE, subblock, "native", MPI_INFO_NULL);
         //reshuffle data
         l=0;
         for(x=0; x<p->geo->lL[1]; x++)
         for(y=0; y<p->geo->lL[2]; y++)
         for(z=0; z<p->geo->lL[3]; z++)
         for(t=0; t<p->geo->lL[0]; t++)
         for(mu=0; mu<4; mu++)
         for(c1=0; c1<3; c1++)
         {
            buffer[l++] = p->D[qcd_LEXIC(t,x,y,z,p->geo->lL)][mu][nu][c1][c2].re;
            buffer[l++] = p->D[qcd_LEXIC(t,x,y,z,p->geo->lL)][mu][nu][c1][c2].im;
         }
         if(!qcd_isBigEndian())      
	   qcd_swap_8(buffer,(size_t)(2*4*3)*(size_t)p->geo->lV);
         
         MPI_File_write_all(mpifid, buffer, 4*3*2*p->geo->lV, MPI_DOUBLE, &status);
         MPI_File_close(&mpifid);
      }
   
   free(buffer);
   MPI_Type_free(&subblock);
//   if(p->geo->myid==0) printf("buffer released\n",c2+nu*3);
   return 0;
}//end qcd_writePropagatorHMCV

int qcd_writePropagator(char *fname, int type, qcd_propagator *p)
{
   if(!p->initialized)
   {
      fprintf(stderr,"Error in qcd_writePropagator! Propagator must be properly initialized.\n");
      return(1);
   }
   switch(type)
   {
      case 0: // CMI format
         //if(p->geo->myid==0) printf("writing CMI type propagator\n");
         return(qcd_writePropagatorCMI(fname, p));
         break;
      case 2:
         return(qcd_writePropagatorHMCV(fname, p));
         break;
      default:
         fprintf(stderr,"process %i: Error in qcd_writePropagator! Unknown type.\n",p->geo->myid);
         return(1);
   }
}//end qcd_writePropagator



int qcd_writeVectorCMI(char *fname, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v)
{
   FILE *fid;
   qcd_uint_4 mu,c1;
   qcd_uint_4 x,y,z,t;
   qcd_uint_8 l;
   int sizes[5],lsizes[5],starts[5];
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   qcd_real_4 *buffer;
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_writeVectorCMI! Vector must be properly initialized.\n");
      return(1);
   }
   
   fid=fopen(fname,"w");
   if(fid==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_writeVectorCMI! Cannot open %s for writing\n",v->geo->myid,fname);
      return(1);
   }
   fclose(fid);
      
   sizes[0]=v->geo->L[1]; 
   sizes[1]=v->geo->L[2];
   sizes[2]=v->geo->L[3];
   sizes[3]=v->geo->L[0];
   sizes[4]=4*3*2;
   lsizes[0]=v->geo->lL[1];
   lsizes[1]=v->geo->lL[2];
   lsizes[2]=v->geo->lL[3];
   lsizes[3]=v->geo->lL[0];
   lsizes[4]=sizes[4];
   starts[0]=v->geo->Pos[1]*lsizes[0];
   starts[1]=v->geo->Pos[2]*lsizes[1];
   starts[2]=v->geo->Pos[3]*lsizes[2];
   starts[3]=v->geo->Pos[0]*lsizes[3];                  
   starts[4]=0;

   buffer = (qcd_real_4*) malloc(v->geo->lV*4*3*2*sizeof(qcd_real_4));
   if(buffer==NULL)
   {
      fprintf(stderr,"process %i: Error in qcd_writeVectorCMI! Out of memory\n",v->geo->myid);
      return(1);
   }
   MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_FLOAT,&subblock);
   MPI_Type_commit(&subblock);
            
   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, 0, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);
   //reshuffle data
   l=0;
   for(x=0; x<v->geo->lL[1]; x++)
   for(y=0; y<v->geo->lL[2]; y++)
   for(z=0; z<v->geo->lL[3]; z++)
   for(t=0; t<v->geo->lL[0]; t++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   {
      buffer[l++] = v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1].re;
      buffer[l++] = v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1].im;
   }
   if(!qcd_isBigEndian())      
     qcd_swap_4(buffer,(size_t)(2*4*3)*(size_t)v->geo->lV);

   MPI_File_write_all(mpifid, buffer, 4*3*2*v->geo->lV, MPI_FLOAT, &status);
   MPI_File_close(&mpifid);         
      
   
   free(buffer);
   MPI_Type_free(&subblock);
   return 0;
}//end qcd_writeVectorCMI



/*writes a double precision ildg formatted vector file*/
int qcd_writeVectorLime(char *fname, int type, qcd_vector *v)
{
   FILE *fid;
   qcd_int_4 i;
   int error_in_header=0;
   LimeWriter *limewriter;
   LimeRecordHeader *limeheader = NULL;
   int ME_flag=0, MB_flag=0, limeStatus;
   qcd_uint_8 message_length;
   MPI_Offset offset;
   MPI_Datatype subblock;  //MPI-type, 5d subarray
   MPI_File mpifid;
   MPI_Status status;
   int sizes[5], lsizes[5], starts[5];
   qcd_uint_2 chunksize,mu,c1;
   char *buffer;
   qcd_uint_4 x,y,z,t;
   char tmp_string[2048];

   if(!v->initialized)
   {
     fprintf(stderr,"Error in %s! Vector must be properly initialized.\n", __func__);
     return(1);
   }   


   //process 0 creates file and writes header   
   if(v->geo->myid == 0)
   {
      fid=fopen(fname,"w");
      if(fid==NULL)
      {
	fprintf(stderr,"process 0: Error in %s! Could not open %s for writing\n",fname, __func__);
	error_in_header=1;
      }else
      {

         limewriter = limeCreateWriter(fid); 
         if(limewriter == (LimeWriter*)NULL) {
	   fprintf(stderr, "Error in %s. LIME error in file %s for writing!\n", fname, __func__);
            error_in_header=1;
         }else
         {            
	   if (type == qcd_PROP_LIME)
	     sprintf(tmp_string, "DiracFermion_Sink");
	   else if (type == qcd_SOURCE_LIME)
	     sprintf(tmp_string, "DiracFermion_Source");

            message_length=(qcd_uint_8) strlen(tmp_string);
            MB_flag=1; ME_flag=1;
	    
	   if (type == qcd_PROP_LIME)
	     limeheader = limeCreateHeader(MB_flag, ME_flag, "propagator-type", message_length);
	   else if (type == qcd_SOURCE_LIME)
	     limeheader = limeCreateHeader(MB_flag, ME_flag, "source-type", message_length);
	   
            if(limeheader == (LimeRecordHeader*)NULL) 
            {
	      fprintf(stderr, "Error in %s. LIME create header error.\n", __func__);
               error_in_header=1;
            }
            limeStatus = limeWriteRecordHeader(limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in %s. LIME write header %d\n", __func__, limeStatus);
              error_in_header=1;
            }
            limeDestroyHeader(limeheader);
            limeStatus = limeWriteRecordData(tmp_string, &message_length, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in %s. LIME write header error %d\n", __func__, limeStatus);
              error_in_header=1;
            }


	    sprintf(tmp_string, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<etmcFormat>\n\t<field>diracFermion</field>\n\t<precision>64</precision>\n\t<flavours>1</flavours>\n\t<lx>%d</lx>\n\t<ly>%d</ly>\n\t<lz>%d</lz>\n\t<lt>%d</lt>\n\t<spin>4</spin>\n\t<colour>3</colour>\n</etmcFormat>", v->geo->L[1], v->geo->L[2], v->geo->L[3], v->geo->L[0]);

            message_length=(qcd_uint_8) strlen(tmp_string);
            MB_flag=1; ME_flag=1;
	    
	    if (type == qcd_PROP_LIME)
	      limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-propagator-format", message_length);
	    else if (type == qcd_SOURCE_LIME)
	      limeheader = limeCreateHeader(MB_flag, ME_flag, "etmc-source-format", message_length);
            if(limeheader == (LimeRecordHeader*)NULL) 
            {
	      fprintf(stderr, "Error in %s. LIME create header error.\n", __func__);
               error_in_header=1;
            }
            limeStatus = limeWriteRecordHeader(limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in %s. LIME write header %d\n", __func__, limeStatus);
              error_in_header=1;
            }
            limeDestroyHeader(limeheader);
            limeStatus = limeWriteRecordData(tmp_string, &message_length, limewriter);
            if(limeStatus < 0 ) 
            {
              fprintf(stderr, "Error in %s. LIME write header error %d\n", __func__, limeStatus);
              error_in_header=1;
            }

            message_length = v->geo->V*4*3*sizeof(qcd_complex_16);
            MB_flag=1; ME_flag=1;
            limeheader = limeCreateHeader(MB_flag, ME_flag, "scidac-binary-data", message_length);
            limeStatus = limeWriteRecordHeader( limeheader, limewriter);
            if(limeStatus < 0 ) 
            {
	      fprintf(stderr, "Error in %s. LIME write header error %d\n", __func__, limeStatus);
	      error_in_header=1;
            }
            limeDestroyHeader( limeheader );
         }
         //make one fake record-write to set offset
         message_length=1;
         limeWriteRecordData(tmp_string, &message_length, limewriter);
         offset = ftell(fid)-1;
         fclose(fid);
      }//fid non null condition
   }//myid==0 condition
   
   MPI_Bcast(&error_in_header,1,MPI_INT,0,MPI_COMM_WORLD);
   if(error_in_header)
     return(1);
   MPI_Bcast(&offset,sizeof(MPI_Offset),MPI_BYTE,0,MPI_COMM_WORLD);
   
   //header written. Now use MPI-I/O to write binary data
   sizes[0]=v->geo->L[0]; 
   sizes[1]=v->geo->L[3];
   sizes[2]=v->geo->L[2];
   sizes[3]=v->geo->L[1];
   sizes[4]=4*3*2;
   lsizes[0]=v->geo->lL[0];
   lsizes[1]=v->geo->lL[3];
   lsizes[2]=v->geo->lL[2];
   lsizes[3]=v->geo->lL[1];
   lsizes[4]=sizes[4];
   starts[0]=v->geo->Pos[0]*lsizes[0];
   starts[1]=v->geo->Pos[3]*lsizes[1];
   starts[2]=v->geo->Pos[2]*lsizes[2];
   starts[3]=v->geo->Pos[1]*lsizes[3];                  
   starts[4]=0;

   MPI_Type_create_subarray(5,sizes,lsizes,starts,MPI_ORDER_C,MPI_DOUBLE,&subblock);
   MPI_Type_commit(&subblock);

   MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_WRONLY, MPI_INFO_NULL, &mpifid);
   MPI_File_set_view(mpifid, offset, MPI_FLOAT, subblock, "native", MPI_INFO_NULL);

   chunksize=4*3*sizeof(qcd_complex_16);
   buffer = (char*) malloc(chunksize*v->geo->lV);
   if(buffer==NULL)
   {
     fprintf(stderr,"process %i: Error in %s! Out of memory\n",v->geo->myid, __func__);
      return(1);
   }
   i=0;
   for(t=0; t<v->geo->lL[0];t++)
   for(z=0; z<v->geo->lL[3];z++)
   for(y=0; y<v->geo->lL[2];y++)
   for(x=0; x<v->geo->lL[1];x++)
     for(mu=0; mu<4; mu++)
       for(c1=0; c1<3; c1++)
	 {
	   ((qcd_real_8 *)buffer)[i] = v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1].re;
	   ((qcd_real_8 *)buffer)[i+1] = v->D[qcd_LEXIC(t,x,y,z,v->geo->lL)][mu][c1].im;
	   i+=2;
	 }
   if(!qcd_isBigEndian())
     qcd_swap_8((double*) buffer,(size_t)(2*4*3)*(size_t)v->geo->lV);
   
   MPI_File_write_all(mpifid, buffer, 4*3*2*v->geo->lV, MPI_DOUBLE, &status);

   free(buffer);
   MPI_File_close(&mpifid);
   MPI_Type_free(&subblock);

   return 0;
}//end qcd_writeVectorLime


int qcd_writeVector(char *fname, int type, qcd_uint_2 nu, qcd_uint_2 c2, qcd_vector *v)
{
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_writeVector! Vector must be properly initialized.\n");
      return(1);
   }
   switch(type)
   {
      case 0: // CMI format
         return(qcd_writeVectorCMI(fname, nu, c2, v));
         break;
      case qcd_PROP_LIME:
      case qcd_SOURCE_LIME:
	return(qcd_writeVectorLime(fname, type, v));
         break;
/*      case 1: // old HMC format
         return(qcd_writeVectorHMC(fname, nu, c2, v));
         break;   */
      default:
         fprintf(stderr,"process %i: Error in qcd_writeVector! Unknown type.\n",v->geo->myid);
         return(1);
   }
}//end qcd_writeVector
