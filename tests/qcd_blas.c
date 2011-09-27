/* qcd_blas.c
 *
 * basic linear algebra on
 * vectors and propagators
 *
 * unoptimized version for all architectures
 *
 * Tomasz Korzec 2008
 ***********************************************/
 
void qcd_zeroVector(qcd_vector *vec)
{
   memset(&(vec->D[0][0][0].re), 0, 4*3*(vec->geo)->lV*sizeof(qcd_complex_16));
}//end qcd_zeroVector

void qcd_setVector(qcd_vector *vec, qcd_complex_16 c)
{
   qcd_uint_8 i;
   for(i=0; i<(vec->geo)->lV; i++)
   {
      vec->D[i][0][0] = c;
      vec->D[i][0][1] = c;
      vec->D[i][0][2] = c;
      vec->D[i][1][0] = c;
      vec->D[i][1][1] = c;
      vec->D[i][1][2] = c;
      vec->D[i][2][0] = c;
      vec->D[i][2][1] = c;
      vec->D[i][2][2] = c;
      vec->D[i][3][0] = c;
      vec->D[i][3][1] = c;
      vec->D[i][3][2] = c;                        
   }
}//end qcd_zeroVector

void qcd_copyVector(qcd_vector *dest, qcd_vector *src)
{
   memcpy(&(dest->D[0][0][0].re),&(src->D[0][0][0].re),(size_t) 4*3*dest->geo->lV*sizeof(qcd_complex_16));
}//end qcd_copyVector

void qcd_copyGaugeField(qcd_gaugeField *dest, qcd_gaugeField *src)
{
   memcpy(&(dest->D[0][0][0][0].re),&(src->D[0][0][0][0].re),(size_t) 4*3*3*dest->geo->lV*sizeof(qcd_complex_16));
}//end qcd_copyGaugeField

void qcd_copyVectorPropagator(qcd_vector *vec, qcd_propagator *prop, qcd_uint_2 nu, qcd_uint_2 c2)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,c1;
   for(i=0; i<vec->geo->lV;i++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
      vec->D[i][mu][c1] = prop->D[i][mu][nu][c1][c2];
}

void qcd_copyVectorPropagator2(qcd_vector *vec, qcd_propagator *prop, qcd_uint_2 mu, qcd_uint_2 c1)
{
   qcd_uint_8 i;
   qcd_uint_2 nu,c2;
   for(i=0; i<vec->geo->lV;i++)
   for(nu=0; nu<4; nu++)
   for(c2=0; c2<3; c2++)
      vec->D[i][nu][c2] = prop->D[i][mu][nu][c1][c2];
}


void qcd_mulVectorC(qcd_vector *vec, qcd_complex_16 c)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a;
   for(i=0; i<vec->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
      vec->D[i][mu][a] = qcd_CMUL(vec->D[i][mu][a],c);
}

/*void qcd_scaleVector(qcd_vector *vec, qcd_real_8 r)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a;
   for(i=0; i<vec->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
      vec->D[i][mu][a] = qcd_CSCALE(vec->D[i][mu][a],r);
}
*/
void qcd_scaleVector(qcd_vector *vec, qcd_real_8 r)
{
   qcd_uint_8 i;
   qcd_real_8 *v;
   v = (qcd_real_8*) &(vec->D[0][0][0].re);
   for(i=0; i<vec->geo->lV*24; i++)
      v[i] *=r;
}



void qcd_scaleVector3d(qcd_vector *vec, qcd_real_8 r, qcd_uint_4 tt)
{
   qcd_uint_8 i;
   qcd_uint_4 t,x,y,z;
   qcd_uint_2 mu,a;
   if(tt<vec->geo->Pos[0]*vec->geo->lL[0] || tt>=(vec->geo->Pos[0]+1)*vec->geo->lL[0])
   {
      return;
   }
   t =tt- vec->geo->Pos[0]*vec->geo->lL[0];
   for(z=0;z<vec->geo->lL[3];z++)
   for(y=0;y<vec->geo->lL[2];y++)
   for(x=0;x<vec->geo->lL[1];x++)
   {
      i=qcd_LEXIC(t,x,y,z,vec->geo->lL);
      for(mu=0; mu<4; mu++)
      for(a=0; a<3; a++)
         vec->D[i][mu][a] = qcd_CSCALE(vec->D[i][mu][a],r);
   }      
}

void qcd_subVector(qcd_vector *difference, qcd_vector *minuend, qcd_vector *subtrahend)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a;
   for(i=0; i<difference->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
      difference->D[i][mu][a] = qcd_CSUB(minuend->D[i][mu][a],subtrahend->D[i][mu][a]);
}

void qcd_subGaugeField(qcd_gaugeField *difference, qcd_gaugeField *minuend, qcd_gaugeField *subtrahend)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a,b;
   for(i=0; i<difference->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
   for(b=0; b<3; b++)
      difference->D[i][mu][a][b] = qcd_CSUB(minuend->D[i][mu][a][b],subtrahend->D[i][mu][a][b]);
}

/*
void qcd_addVector(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a;
   for(i=0; i<sum->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
      sum->D[i][mu][a] = qcd_CADD(summand1->D[i][mu][a],summand2->D[i][mu][a]);
}
*/
void qcd_addVector(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2)
{
   qcd_uint_8 i;
   qcd_real_8 *su, *s1, *s2;
   su = (qcd_real_8*) &(sum->D[0][0][0].re);
   s1 = (qcd_real_8*) &(summand1->D[0][0][0].re);
   s2 = (qcd_real_8*) &(summand2->D[0][0][0].re);      
   for(i=0; i<sum->geo->lV*24; i++)
      su[i] = s1[i]+s2[i];
}


void qcd_addVector3d(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2, qcd_uint_4 tt)
{
   qcd_uint_8 i;
   qcd_uint_4 t,x,y,z;
   qcd_uint_2 mu,a;
   if(tt<sum->geo->Pos[0]*sum->geo->lL[0] || tt>=(sum->geo->Pos[0]+1)*sum->geo->lL[0])
   {
      return; //this process doesn't contain the time-slice
   }
   t = tt-sum->geo->Pos[0]*sum->geo->lL[0];
   for(z=0;z<sum->geo->lL[3];z++)
   for(y=0;y<sum->geo->lL[2];y++)
   for(x=0;x<sum->geo->lL[1];x++)
   {
      i=qcd_LEXIC(t,x,y,z,sum->geo->lL);
      for(mu=0; mu<4; mu++)
      for(a=0; a<3; a++)
         sum->D[i][mu][a] = qcd_CADD(summand1->D[i][mu][a],summand2->D[i][mu][a]);
   }      
}

double qcd_normVector(qcd_vector *vec)
{
   qcd_real_8 tmp=0, result=0;
   qcd_uint_8 i;
   qcd_uint_2 mu,a;
      
   for(i=0; i<vec->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
      tmp += qcd_NORM(vec->D[i][mu][a]);
      
   MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return((double) result);   
}

/* sum of norms of all links. norm(link) = | trace(link) 
 */
double qcd_normGaugeField(qcd_gaugeField *u)
{
   qcd_real_8 tmp=0, result=0;
   qcd_uint_8 i;
   qcd_uint_2 mu,a;
   qcd_complex_16 tr;
      
   for(i=0; i<u->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   {
      tr = (qcd_complex_16) {0,0};
      for(a=0; a<3; a++)
         tr = qcd_CADD(tr,u->D[i][mu][a][a]);
   tmp += qcd_NORM(tr);    
   }
   MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return((double) result);   
}
