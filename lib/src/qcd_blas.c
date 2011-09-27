/* qcd_blas.c
 *
 * basic linear algebra on
 * vectors and propagators
 *
 * unoptimized version for all architectures
 *
 * Tomasz Korzec 2008
 ***********************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>


void qcd_zero3x3(qcd_complex_16 a[3][3])
{
   memset(&(a[0][0].re),0,9*sizeof(qcd_complex_16));
}//end qcd_zero3x3

void qcd_unit3x3(qcd_complex_16 a[3][3])
{
   memset(&(a[0][0].re),0,9*sizeof(qcd_complex_16));
   a[0][0].re=1.0;
   a[1][1].re=1.0;
   a[2][2].re=1.0;
   return;
}//end qcd_unit3x3

void qcd_minusUnit3x3(qcd_complex_16 a[3][3])
{
   memset(&(a[0][0].re),0,9*sizeof(qcd_complex_16));
   a[0][0].re=-1.0;
   a[1][1].re=-1.0;
   a[2][2].re=-1.0;
   return;
}//end qcd_minusUnit3x3

void qcd_sub3x3(qcd_complex_16 difference[3][3], qcd_complex_16 minuend[3][3], qcd_complex_16 subtrahend[3][3])
{ 
   qcd_uint_2 c1,c2;
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
      difference[c1][c2] = qcd_CSUB(minuend[c1][c2],subtrahend[c1][c2]);
   return;   
}//end qcd_sub3x3

void qcd_add3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3])
{ 
   qcd_uint_2 c1,c2;
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
      sum[c1][c2] = qcd_CADD(summand1[c1][c2],summand2[c1][c2]);
   return;   
}//end qcd_add3x3

void qcd_addAdjoint3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3])
{ 
   qcd_uint_2 c1,c2;
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
      sum[c1][c2] = qcd_CADD(summand1[c1][c2],qcd_CONJ(summand2[c2][c1]));
   return;   
}//end qcd_addAdjoint3x3

void qcd_mul3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3])
{ 
   qcd_uint_2 c1,c2,c3;
   memset(&(product[0][0].re),0,9*sizeof(qcd_complex_16));
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
   for(c3=0;c3<3;c3++)
   {
      product[c1][c2] = qcd_CADD(product[c1][c2], qcd_CMUL(factor1[c1][c3],factor2[c3][c2]));
   }   
   return;   
}//end qcd_mul3x3

void qcd_mulAdjoint3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3])
{ 
   qcd_uint_2 c1,c2,c3;
   memset(&(product[0][0].re),0,9*sizeof(qcd_complex_16));
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
   for(c3=0;c3<3;c3++)
   {
      product[c1][c2] = qcd_CADD(product[c1][c2], qcd_CMUL(factor1[c1][c3],qcd_CONJ(factor2[c2][c3])));
   }   
   return;   
}//end qcd_mulAdjoint3x3

void qcd_scale3x3(qcd_complex_16 a[3][3], qcd_real_8 r)
{ 
   qcd_uint_2 c1,c2;
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
      a[c1][c2] = qcd_CSCALE(a[c1][c2], r);
   return;   
}//end qcd_scale3x3

void qcd_cScale3x3(qcd_complex_16 a[3][3], qcd_complex_16 c)
{ 
   qcd_uint_2 c1,c2;
   for(c1=0;c1<3;c1++)
   for(c2=0;c2<3;c2++)
      a[c1][c2] = qcd_CMUL(a[c1][c2], c);
   return;   
}//end qcd_cScale3x3

qcd_complex_16 qcd_trace3x3(qcd_complex_16 a[3][3])
{
   qcd_complex_16 res;
   res.re = a[0][0].re + a[1][1].re + a[2][2].re;
   res.im = a[0][0].im + a[1][1].im + a[2][2].im;
   return(res);
}//end qcd_trace3x3

void qcd_copy3x3(qcd_complex_16 dest[3][3], qcd_complex_16 src[3][3])
{
   memcpy(&(dest[0][0].re),&(src[0][0].re),(size_t) 9*sizeof(qcd_complex_16));
}//end qcd_copy3x3

void qcd_dagger3x3(qcd_complex_16 a[3][3])
{
   qcd_complex_16 tmp[3][3];
   qcd_uint_2 c1,c2;
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
      tmp[c1][c2] = qcd_CONJ(a[c2][c1]);
   memcpy(&(a[0][0].re),&(tmp[0][0].re),(size_t) 9*sizeof(qcd_complex_16));
   return;   
}//end qcd_dagger3x3
 
void qcd_zeroVector(qcd_vector *vec)
{
   memset(&(vec->D[0][0][0].re), 0, 4*3*(vec->geo)->lV*sizeof(qcd_complex_16));
}//end qcd_zeroVector

void qcd_zeroPropagator(qcd_propagator *prop)
{
   memset(&(prop->D[0][0][0][0][0].re), 0, 4*4*3*3*(prop->geo)->lV*sizeof(qcd_complex_16));
}//end qcd_zeroPropagator

void qcd_zeroGaugeField(qcd_gaugeField *u)
{
   memset(&(u->D[0][0][0][0].re), 0, 4*3*3*(u->geo)->lV*sizeof(qcd_complex_16));
}//end qcd_zeroGaugeField

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

void qcd_copyPropagatorVector(qcd_propagator *prop, qcd_vector *vec, qcd_uint_2 nu, qcd_uint_2 c2)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,c1;
   for(i=0; i<vec->geo->lV;i++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
      prop->D[i][mu][nu][c1][c2] = vec->D[i][mu][c1];
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


void qcd_copyPropagatorVector2(qcd_propagator *prop, qcd_vector *vec, qcd_uint_2 mu, qcd_uint_2 c1)
{
   qcd_uint_8 i;
   qcd_uint_2 nu,c2;
   for(i=0; i<vec->geo->lV;i++)
   for(nu=0; nu<4; nu++)
   for(c2=0; c2<3; c2++)
      prop->D[i][mu][nu][c1][c2] = vec->D[i][nu][c2];
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


void qcd_mulPropagatorC(qcd_propagator *prop, qcd_complex_16 c)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,nu,a,b;
   for(i=0; i<prop->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   for(a=0; a<3; a++)
   for(b=0; b<3; b++)
      prop->D[i][mu][nu][a][b] = qcd_CMUL(prop->D[i][mu][nu][a][b],c);
}


void qcd_mulPropagatorC3d(qcd_propagator *prop, qcd_complex_16 c, qcd_uint_4 t)
{
   qcd_uint_8 i,v;
   qcd_uint_2 mu,nu,a,b;
   for(i=0; i<prop->geo->lV3; i++)
   {
      v = i*prop->geo->lL[0] + t; // works only with current lexicographical conventions
      for(mu=0; mu<4; mu++)
      for(nu=0; nu<4; nu++)
      for(a=0; a<3; a++)
      for(b=0; b<3; b++)
         prop->D[v][mu][nu][a][b] = qcd_CMUL(prop->D[v][mu][nu][a][b],c);
   }      
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

/*
void qcd_subVector(qcd_vector *difference, qcd_vector *minuend, qcd_vector *subtrahend)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a;
   for(i=0; i<difference->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
      difference->D[i][mu][a] = qcd_CSUB(minuend->D[i][mu][a],subtrahend->D[i][mu][a]);
}
*/

void qcd_subVector(qcd_vector *difference, qcd_vector *minuend, qcd_vector *subtrahend)
{
   qcd_uint_8 i;
   qcd_real_8 *su, *s1, *s2;
   su = (qcd_real_8*) &(difference->D[0][0][0].re);
   s1 = (qcd_real_8*) &(minuend->D[0][0][0].re);
   s2 = (qcd_real_8*) &(subtrahend->D[0][0][0].re);      
   for(i=0; i<difference->geo->lV*24; i++)
      su[i] = s1[i]-s2[i];
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

void qcd_addGaugeField(qcd_gaugeField *sum, qcd_gaugeField *summand1, qcd_gaugeField *summand2)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a,b;
   for(i=0; i<sum->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
   for(b=0; b<3; b++)
      sum->D[i][mu][a][b] = qcd_CADD(summand1->D[i][mu][a][b],summand2->D[i][mu][a][b]);
}

void qcd_scaleGaugeField(qcd_gaugeField *u, qcd_real_8 alpha)
{
   qcd_uint_8 i;
   qcd_uint_2 mu,a,b;
   for(i=0; i<u->geo->lV; i++)
   for(mu=0; mu<4; mu++)
   for(a=0; a<3; a++)
   for(b=0; b<3; b++)
      u->D[i][mu][a][b] = qcd_CSCALE(u->D[i][mu][a][b],alpha);
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

void qcd_axpyVector(qcd_vector *axpy, qcd_complex_16 a, qcd_vector *x, qcd_vector *y)
{
   /* axpy = a*x + y, a complex, x,y vectors */
   qcd_uint_8 i;
   qcd_real_8 *su, *s1, *s2;
   qcd_real_8 ra,ia;

   ra = a.re;
   ia = a.im;
   su = (qcd_real_8*) &(axpy->D[0][0][0].re);
   s1 = (qcd_real_8*) &(x->D[0][0][0].re);
   s2 = (qcd_real_8*) &(y->D[0][0][0].re);      
   for(i=0; i<x->geo->lV*24; i++)
   {
      su[i] = s2[i]   + s1[i]*ra - s1[i+1]*ia;
      i++;
      su[i] = s2[i]   + s1[i]*ra + s1[i-1]*ia;
   }
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


qcd_complex_16 qcd_mulAdjointVectorVector(qcd_vector *factor1, qcd_vector *factor2)
{
   qcd_uint_4 i;
   qcd_complex_16 *ptr1, *ptr2;
   qcd_complex_16 product= {0,0};
   qcd_complex_16 result= {0,0};
   
   ptr1 = &(factor1->D[0][0][0]);
   ptr2 = &(factor2->D[0][0][0]);
   
   for(i=0; i<factor1->geo->lV*12; i++)
   {
      product.re += ptr1->re * ptr2->re + ptr1->im * ptr2->im;
      product.im += ptr1->re * ptr2->im - ptr1->im * ptr2->re;
      ptr1++;
      ptr2++;
   }
   MPI_Allreduce(&(product.re), &(result.re), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   return(result);
}


qcd_real_8 qcd_normVector(qcd_vector *vec)
{
   qcd_real_8 tmp=0, result=0;
   qcd_uint_8 i;
   qcd_real_8 *v;
      
   v = &(vec->D[0][0][0].re);
   for(i=0; i<vec->geo->lV*24; i++)
      tmp += v[i]*v[i];
   
   MPI_Allreduce(&tmp, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   result = sqrt(result);
   return(result);   
}

/* sum of norms of all links. qcd_NORM(link) = | trace(link) 
 */
qcd_real_8 qcd_normGaugeField(qcd_gaugeField *u)
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
   return(result);   
}



void qcd_conjPropagator(qcd_propagator *prop)
{
   qcd_uint_2 mu,nu,c1,c2;
   qcd_uint_8 l;
   for(l=0; l<prop->geo->lV; l++)
   for(mu=0; mu<4; mu++)
   for(nu=0; nu<4; nu++)
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      prop->D[l][mu][nu][c1][c2] = qcd_CONJ(prop->D[l][mu][nu][c1][c2]);
   }
}



/*
  Projects the gauge_Field to SU(3). Written by Giannis Koutsou. Adapted to qcd_lib format by T.K.
  If M is an arbitary 3x3 matrix it can be factorised as the product of a hermitian and a unitary
  matrix, M = U H. We will calculate U as the projection of M to SU(3). Now M^\dagger M = H^2. By
  going to diagonal basis for H^2 we calculate H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^\dagger where v
  is the 3x3 matrix whoes columns are the eigenvectors of H^2 and L is the diagonal representation
  of H^2. Thus: U = M H^-1
  [1] Phys. Lett. B307, 375-382, (1993)
*/
void qcd_projectSU3(qcd_gaugeField *gf)
{
   qcd_uint_4 iv,mu,c1,c2;
   qcd_complex_16 H[3][3],detM,U[3][3],M[3][3];
   qcd_complex_16 b,D,v[3][3],vr[3][3];
   qcd_real_8 a,ThirdRoot_18,ThirdRoot_12,ThirdRoot_2_3;
   qcd_real_8 trace,e[3],de,cor[3];
   qcd_complex_16 w,ThirdRootOne[2];
   qcd_real_8  sum;
   qcd_real_8  norm,phase;

   /*
     Constants Used:
   */
   ThirdRootOne[0].re= 1;
   ThirdRootOne[0].im= sqrt(3.);

   ThirdRootOne[1].re= 1;
   ThirdRootOne[1].im=-sqrt(3.);

   ThirdRoot_12 =pow(12.,1./3.);
   ThirdRoot_18 =pow(18.,1./3.);
   ThirdRoot_2_3=pow((2./3.),1./3.);

   for(iv=0; iv<gf->geo->lV; iv++)
   for(mu=0; mu<4; mu++) 
   {
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
         M[c1][c2] = gf->D[iv][mu][c1][c2];

      
      detM = qcd_CADD(qcd_CADD( qcd_CMUL(M[0][0],qcd_CMUL(M[1][1],M[2][2])),
                                qcd_CMUL(M[0][1],qcd_CMUL(M[1][2],M[2][0]))),
                                qcd_CMUL(M[0][2],qcd_CMUL(M[1][0],M[2][1])));
                                
      detM = qcd_CSUB(detM,
             qcd_CADD(qcd_CADD( qcd_CMUL(M[0][2],qcd_CMUL(M[1][1],M[2][0])),
                                qcd_CMUL(M[0][0],qcd_CMUL(M[1][2],M[2][1]))),
                                qcd_CMUL(M[0][1],qcd_CMUL(M[1][0],M[2][2]))));

      phase = qcd_ARG(detM)/3.;

      H[0][0].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][0]);
      H[0][1].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][1]);
      H[0][2].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][2]);
      H[0][0].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][0]);
      H[0][1].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][1]);
      H[0][2].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][2]);

      H[1][0].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][0]);
      H[1][1].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][1]);
      H[1][2].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][2]);
      H[1][0].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][0]);
      H[1][1].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][1]);
      H[1][2].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][2]);

      H[2][0].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][0]);
      H[2][1].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][1]);
      H[2][2].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][2]);
      H[2][0].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][0]);
      H[2][1].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][1]);
      H[2][2].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][2]);

      /*
        Assure Hermiticity:
      */
      H[0][1].re = (H[0][1].re + H[1][0].re)/2.;
      H[0][1].im = (H[0][1].im - H[1][0].im)/2.;

      H[1][0] =  qcd_CONJ(H[0][1]);

      H[0][2].re = (H[0][2].re + H[2][0].re)/2.;
      H[0][2].im = (H[0][2].im - H[2][0].im)/2.;

      H[2][0] =  qcd_CONJ(H[0][2]);

      H[1][2].re = (H[1][2].re + H[2][1].re)/2.;
      H[1][2].im = (H[1][2].im - H[2][1].im)/2.;

      H[2][1] =  qcd_CONJ(H[1][2]);
      /*
        If H^2 is alread diagonal skip diagonalization and
        calculate U directly
      */
      sum=qcd_NORM(H[0][1])+qcd_NORM(H[0][2])+qcd_NORM(H[1][2]);

      if(sum<=1e-08)
      {
         e[0]=1./sqrt(H[0][0].re);
         e[1]=1./sqrt(H[1][1].re);
         e[2]=1./sqrt(H[2][2].re);

         U[0][0] = (qcd_complex_16) { M[0][0].re*e[0], M[0][0].im*e[0] };
         U[0][1] = (qcd_complex_16) { M[0][1].re*e[0], M[0][1].im*e[0] };
         U[0][2] = (qcd_complex_16) { M[0][2].re*e[0], M[0][2].im*e[0] };

         U[1][0] = (qcd_complex_16) { M[1][0].re*e[1], M[1][0].im*e[1] };
         U[1][1] = (qcd_complex_16) { M[1][1].re*e[1], M[1][1].im*e[1] };
         U[1][2] = (qcd_complex_16) { M[1][2].re*e[1], M[1][2].im*e[1] };

         U[2][0] = (qcd_complex_16) { M[2][0].re*e[2], M[2][0].im*e[2] };
         U[2][1] = (qcd_complex_16) { M[2][1].re*e[2], M[2][1].im*e[2] };
         U[2][2] = (qcd_complex_16) { M[2][2].re*e[2], M[2][2].im*e[2] };

         for(c1=0; c1<3; c1++)
         for(c2=0; c2<3; c2++)
            gf->D[iv][mu][c1][c2] = U[c1][c2]; 
      }
      else
      {
         /*
           Make traceless to eliminate second order term in eigenvalue equation,
           i.e. eig^3 + a eig + b = 0, when H is traceless.
         */
         trace=(H[0][0].re+H[1][1].re+H[2][2].re)/3.;

         H[0][0].re-=trace;
         H[1][1].re-=trace;
         H[2][2].re-=trace;


         /*
           Solve for eigenvalues:
           e^3 - e (H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2) - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 + H11*|H23|^2 - H13*H12^* *H23^*,
           e^3 + a*e + b = 0,

           a = -(H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2)
           b = - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 - H22*|H23|^2 - H33*|H23|^2 - H13*H12^* *H23^*

           D=(-9b + sqrt(12a^3 + 81b^2))^(1/3)

           e = D/(18^(1/3)) - ((2/3)^(1/3))/D
           e = (1 + I sqrt(3))a / (D 12^(1/3)) - (1 - I sqrt(3)) D / (2 18^(1/3))
           e = (1 - I sqrt(3))a / (D 12^(1/3)) - (1 + I sqrt(3)) D / (2 18^(1/3))
         */
         a = -(H[2][2].re*H[2][2].re - H[0][0].re*H[1][1].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + qcd_CMULR(H[0][2],qcd_CONJ(H[0][2])) + qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])));

         b.re  = - H[0][0].re*H[1][1].re*H[2][2].re + H[2][2].re*qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULR(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULR(H[0][2],qcd_CONJ(H[0][2]));
         b.im  =                                      H[2][2].re*qcd_CMULI(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULI(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULI(H[0][2],qcd_CONJ(H[0][2]));

         b.re +=   H[0][0].re*qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULR(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));
         b.im +=   H[0][0].re*qcd_CMULI(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULI(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));

         w.re=qcd_CPOWR(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);
         w.im=qcd_CPOWI(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);

         D=qcd_CPOW(((qcd_complex_16){-9.*b.re + w.re,-9.*b.im + w.im}),1./3.); 

         e[0] = D.re/(ThirdRoot_18) - qcd_CDEVR(((qcd_complex_16){a*ThirdRoot_2_3,0}),D);
         e[1] = a*qcd_CDEVR(ThirdRootOne[0],((qcd_complex_16){D.re*ThirdRoot_12,D.im*ThirdRoot_12})) - qcd_CMULR(ThirdRootOne[1],D)/(ThirdRoot_18*2.);
         e[2] = -e[0]-e[1];

         e[0]+= trace;
         e[1]+= trace;
         e[2]+= trace;

         H[0][0].re+=trace;
         H[1][1].re+=trace;
         H[2][2].re+=trace;

         /*
           Eigenvectors:
           v[0] = -(e H31 - H31 H22 + H21 H32) v[2] / Denom
           v[1] = -(H31 H12 - e H32 - H11 H32) v[2] / Denom
           v[2] =  (-e^2) + e H11 + |H12|^2 + e H22 - H11 H22
         */

         v[0][0].re = -(e[0]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[0][0].im = -(e[0]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[0][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[0]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[0][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[0]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[0][2].re =-e[0]*e[0] + e[0]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[0]*H[1][1].re - H[0][0].re*H[1][1].re;
         v[0][2].im = 0.;

         v[1][0].re = -(e[1]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[1][0].im = -(e[1]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[1][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[1]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[1][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[1]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[1][2].re =-e[1]*e[1] + e[1]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[1]*H[1][1].re - H[0][0].re*H[1][1].re;;
         v[1][2].im = 0.;

         /*
           Assure eigenvectors orthonormality:
           norm =  inner product v1.v1
           w    = (inner product v1.v2)/norm
           v2   = w*v1
         */

         norm  = qcd_CMULR(v[0][0],qcd_CONJ(v[0][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[0][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[0][2]));
         w.re  = qcd_CMULR(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[1][2]));
         w.im  = qcd_CMULI(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULI(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULI(v[0][2],qcd_CONJ(v[1][2]));
         w.re /= norm;
         w.im /= norm;

         v[1][0].re-= qcd_CMULR(w,v[0][0]);
         v[1][0].im-= qcd_CMULI(w,v[0][0]);

         v[1][1].re-= qcd_CMULR(w,v[0][1]);
         v[1][1].im-= qcd_CMULI(w,v[0][1]);

         v[1][2].re-= qcd_CMULR(w,v[0][2]);
         v[1][2].im-= qcd_CMULI(w,v[0][2]);

         norm=1./sqrt(norm);

         /*
           Normalize first and second eigenvector:
         */

         v[0][0].re*= norm;
         v[0][0].im*= norm;

         v[0][1].re*= norm;
         v[0][1].im*= norm;

         v[0][2].re*= norm;
         v[0][2].im*= norm;


         norm = qcd_CMULR(v[1][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[1][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[1][2],qcd_CONJ(v[1][2]));

         norm=1./sqrt(norm);

         v[1][0].re*= norm;
         v[1][0].im*= norm;

         v[1][1].re*= norm;
         v[1][1].im*= norm;

         v[1][2].re*= norm;
         v[1][2].im*= norm;

         /*
           v3 = v1 x v2
         */


         v[2][0].re =  qcd_CMULR(v[0][1],v[1][2]) - qcd_CMULR(v[0][2],v[1][1]);
         v[2][0].im = -qcd_CMULI(v[0][1],v[1][2]) + qcd_CMULI(v[0][2],v[1][1]);

         v[2][1].re = -qcd_CMULR(v[0][0],v[1][2]) + qcd_CMULR(v[0][2],v[1][0]);
         v[2][1].im = +qcd_CMULI(v[0][0],v[1][2]) - qcd_CMULI(v[0][2],v[1][0]);

         v[2][2].re =  qcd_CMULR(v[0][0],v[1][1]) - qcd_CMULR(v[0][1],v[1][0]);
         v[2][2].im = -qcd_CMULI(v[0][0],v[1][1]) + qcd_CMULI(v[0][1],v[1][0]);

         de     =               e[0]*e[1] +   e[1]*e[2] +   e[2]*e[0];
         /*
         cor[0] = tan(phase) * (e[0]*e[1] - 2*e[1]*e[2] +   e[2]*e[0])/de;
         cor[1] = tan(phase) * (e[0]*e[1] +   e[1]*e[2] - 2*e[2]*e[0])/de;
         cor[2] = - cor[0] - cor[1];
         */
         //to be compatible with Grenoble & Paris, don't apply corrections
         cor[0]=0;
         cor[1]=0;
         cor[2]=0;

         de = 1./sqrt(e[0]);
         b.re = de*cos(phase-cor[0]);
         b.im =-de*sin(phase-cor[0]);
         vr[0][0] = qcd_CMUL(b,v[0][0]);
         vr[0][1] = qcd_CMUL(b,v[0][1]);
         vr[0][2] = qcd_CMUL(b,v[0][2]);

         de = 1./sqrt(e[1]);
         b.re = de*cos(phase-cor[1]);
         b.im =-de*sin(phase-cor[1]);

         vr[1][0] = qcd_CMUL(b,v[1][0]);
         vr[1][1] = qcd_CMUL(b,v[1][1]);
         vr[1][2] = qcd_CMUL(b,v[1][2]);

         de = 1./sqrt(e[2]);
         b.re = de*cos(phase-cor[2]);
         b.im =-de*sin(phase-cor[2]);

         vr[2][0] = qcd_CMUL(b,v[2][0]);
         vr[2][1] = qcd_CMUL(b,v[2][1]);
         vr[2][2] = qcd_CMUL(b,v[2][2]);


         H[0][0].re= qcd_CMULR(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].re= qcd_CMULR(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].re= qcd_CMULR(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[2][2])) ;

         H[0][0].im= qcd_CMULI(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].im= qcd_CMULI(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].im= qcd_CMULI(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[2][2])) ;


         H[1][0].re= qcd_CMULR(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].re= qcd_CMULR(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].re= qcd_CMULR(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[2][2])) ;

         H[1][0].im= qcd_CMULI(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].im= qcd_CMULI(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].im= qcd_CMULI(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[2][2])) ;


         H[2][0].re= qcd_CMULR(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].re= qcd_CMULR(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].re= qcd_CMULR(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[2][2])) ;

         H[2][0].im= qcd_CMULI(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].im= qcd_CMULI(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].im= qcd_CMULI(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[2][2])) ;

         U[0][0].re= qcd_CMULR(H[0][0],vr[0][0]) +  qcd_CMULR(H[0][1],vr[1][0])  +  qcd_CMULR(H[0][2],vr[2][0]) ;
         U[0][1].re= qcd_CMULR(H[0][0],vr[0][1]) +  qcd_CMULR(H[0][1],vr[1][1])  +  qcd_CMULR(H[0][2],vr[2][1]) ;
         U[0][2].re= qcd_CMULR(H[0][0],vr[0][2]) +  qcd_CMULR(H[0][1],vr[1][2])  +  qcd_CMULR(H[0][2],vr[2][2]) ;

         U[0][0].im= qcd_CMULI(H[0][0],vr[0][0]) +  qcd_CMULI(H[0][1],vr[1][0])  +  qcd_CMULI(H[0][2],vr[2][0]) ;
         U[0][1].im= qcd_CMULI(H[0][0],vr[0][1]) +  qcd_CMULI(H[0][1],vr[1][1])  +  qcd_CMULI(H[0][2],vr[2][1]) ;
         U[0][2].im= qcd_CMULI(H[0][0],vr[0][2]) +  qcd_CMULI(H[0][1],vr[1][2])  +  qcd_CMULI(H[0][2],vr[2][2]) ;


         U[1][0].re= qcd_CMULR(H[1][0],vr[0][0]) +  qcd_CMULR(H[1][1],vr[1][0])  +  qcd_CMULR(H[1][2],vr[2][0]) ;
         U[1][1].re= qcd_CMULR(H[1][0],vr[0][1]) +  qcd_CMULR(H[1][1],vr[1][1])  +  qcd_CMULR(H[1][2],vr[2][1]) ;
         U[1][2].re= qcd_CMULR(H[1][0],vr[0][2]) +  qcd_CMULR(H[1][1],vr[1][2])  +  qcd_CMULR(H[1][2],vr[2][2]) ;

         U[1][0].im= qcd_CMULI(H[1][0],vr[0][0]) +  qcd_CMULI(H[1][1],vr[1][0])  +  qcd_CMULI(H[1][2],vr[2][0]) ;
         U[1][1].im= qcd_CMULI(H[1][0],vr[0][1]) +  qcd_CMULI(H[1][1],vr[1][1])  +  qcd_CMULI(H[1][2],vr[2][1]) ;
         U[1][2].im= qcd_CMULI(H[1][0],vr[0][2]) +  qcd_CMULI(H[1][1],vr[1][2])  +  qcd_CMULI(H[1][2],vr[2][2]) ;


         U[2][0].re= qcd_CMULR(H[2][0],vr[0][0]) +  qcd_CMULR(H[2][1],vr[1][0])  +  qcd_CMULR(H[2][2],vr[2][0]) ;
         U[2][1].re= qcd_CMULR(H[2][0],vr[0][1]) +  qcd_CMULR(H[2][1],vr[1][1])  +  qcd_CMULR(H[2][2],vr[2][1]) ;
         U[2][2].re= qcd_CMULR(H[2][0],vr[0][2]) +  qcd_CMULR(H[2][1],vr[1][2])  +  qcd_CMULR(H[2][2],vr[2][2]) ;

         U[2][0].im= qcd_CMULI(H[2][0],vr[0][0]) +  qcd_CMULI(H[2][1],vr[1][0])  +  qcd_CMULI(H[2][2],vr[2][0]) ;
         U[2][1].im= qcd_CMULI(H[2][0],vr[0][1]) +  qcd_CMULI(H[2][1],vr[1][1])  +  qcd_CMULI(H[2][2],vr[2][1]) ;
         U[2][2].im= qcd_CMULI(H[2][0],vr[0][2]) +  qcd_CMULI(H[2][1],vr[1][2])  +  qcd_CMULI(H[2][2],vr[2][2]) ;


         /*
           w    = inner product: col1.col2
           norm = inner product: col1.col1
         */

         norm  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][0])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][0])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][0]));
         w.re  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][1]));
         w.im  = qcd_CMULI(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULI(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULI(U[2][0],qcd_CONJ(U[2][1]));
         w.re /= norm;
         w.im /= norm;


         U[0][1].re-=qcd_CMULR(w,U[0][0]);
         U[0][1].im-=qcd_CMULI(w,U[0][0]);

         U[1][1].re-=qcd_CMULR(w,U[1][0]);
         U[1][1].im-=qcd_CMULI(w,U[1][0]);

         U[2][1].re-=qcd_CMULR(w,U[2][0]);
         U[2][1].im-=qcd_CMULI(w,U[2][0]);

         norm = 1./sqrt(norm);

         U[0][0].re*= norm;
         U[0][0].im*= norm;
         U[1][0].re*= norm;
         U[1][0].im*= norm;
         U[2][0].re*= norm;
         U[2][0].im*= norm;

         norm = qcd_CMULR(U[0][1],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][1],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][1],qcd_CONJ(U[2][1]));
         norm = 1./sqrt(norm);

         U[0][1].re*= norm;
         U[0][1].im*= norm;
         U[1][1].re*= norm;
         U[1][1].im*= norm;
         U[2][1].re*= norm;
         U[2][1].im*= norm;

         /*
           col3 = col1 x col2
         */
         U[0][2].re =  qcd_CMULR(U[1][0],U[2][1]) - qcd_CMULR(U[2][0],U[1][1]);
         U[0][2].im = -qcd_CMULI(U[1][0],U[2][1]) + qcd_CMULI(U[2][0],U[1][1]);

         U[1][2].re = -qcd_CMULR(U[0][0],U[2][1]) + qcd_CMULR(U[2][0],U[0][1]);
         U[1][2].im =  qcd_CMULI(U[0][0],U[2][1]) - qcd_CMULI(U[2][0],U[0][1]);

         U[2][2].re =  qcd_CMULR(U[0][0],U[1][1]) - qcd_CMULR(U[1][0],U[0][1]);
         U[2][2].im = -qcd_CMULI(U[0][0],U[1][1]) + qcd_CMULI(U[1][0],U[0][1]);

         for(c1=0; c1<3; c1++)
         for(c2=0; c2<3; c2++)
            gf->D[iv][mu][c1][c2] = U[c1][c2]; 
      }
   }//end volume-mu-loop
   return;
}//end qcd_projectSU3

/*
  Projects the gauge_Field to SU(3). Written by Giannis Koutsou. Adapted to qcd_lib format by T.K.
  If M is an arbitary 3x3 matrix it can be factorised as the product of a hermitian and a unitary
  matrix, M = U H. We will calculate U as the projection of M to SU(3). Now M^\dagger M = H^2. By
  going to diagonal basis for H^2 we calculate H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^\dagger where v
  is the 3x3 matrix whoes columns are the eigenvectors of H^2 and L is the diagonal representation
  of H^2. Thus: U = M H^-1
  [1] Phys. Lett. B307, 375-382, (1993)
  version that projects spatial links only
*/
void qcd_projectSU33d(qcd_gaugeField *gf)
{
   qcd_uint_4 iv,mu,c1,c2;
   qcd_complex_16 H[3][3],detM,U[3][3],M[3][3];
   qcd_complex_16 b,D,v[3][3],vr[3][3];
   qcd_real_8 a,ThirdRoot_18,ThirdRoot_12,ThirdRoot_2_3;
   qcd_real_8 trace,e[3],de,cor[3];
   qcd_complex_16 w,ThirdRootOne[2];
   qcd_real_8  sum;
   qcd_real_8  norm,phase;

   /*
     Constants Used:
   */
   ThirdRootOne[0].re= 1;
   ThirdRootOne[0].im= sqrt(3.);

   ThirdRootOne[1].re= 1;
   ThirdRootOne[1].im=-sqrt(3.);

   ThirdRoot_12 =pow(12.,1./3.);
   ThirdRoot_18 =pow(18.,1./3.);
   ThirdRoot_2_3=pow((2./3.),1./3.);

   for(iv=0; iv<gf->geo->lV; iv++)
   for(mu=1; mu<4; mu++) 
   {
      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
         M[c1][c2] = gf->D[iv][mu][c1][c2];

      
      detM = qcd_CADD(qcd_CADD( qcd_CMUL(M[0][0],qcd_CMUL(M[1][1],M[2][2])),
                                qcd_CMUL(M[0][1],qcd_CMUL(M[1][2],M[2][0]))),
                                qcd_CMUL(M[0][2],qcd_CMUL(M[1][0],M[2][1])));
                                
      detM = qcd_CSUB(detM,
             qcd_CADD(qcd_CADD( qcd_CMUL(M[0][2],qcd_CMUL(M[1][1],M[2][0])),
                                qcd_CMUL(M[0][0],qcd_CMUL(M[1][2],M[2][1]))),
                                qcd_CMUL(M[0][1],qcd_CMUL(M[1][0],M[2][2]))));

      phase = qcd_ARG(detM)/3.;

      H[0][0].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][0]);
      H[0][1].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][1]);
      H[0][2].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][2]);
      H[0][0].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][0]);
      H[0][1].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][1]);
      H[0][2].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][2]);

      H[1][0].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][0]);
      H[1][1].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][1]);
      H[1][2].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][2]);
      H[1][0].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][0]);
      H[1][1].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][1]);
      H[1][2].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][2]);

      H[2][0].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][0]);
      H[2][1].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][1]);
      H[2][2].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][2]);
      H[2][0].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][0]);
      H[2][1].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][1]);
      H[2][2].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][2]);

      /*
        Assure Hermiticity:
      */
      H[0][1].re = (H[0][1].re + H[1][0].re)/2.;
      H[0][1].im = (H[0][1].im - H[1][0].im)/2.;

      H[1][0] =  qcd_CONJ(H[0][1]);

      H[0][2].re = (H[0][2].re + H[2][0].re)/2.;
      H[0][2].im = (H[0][2].im - H[2][0].im)/2.;

      H[2][0] =  qcd_CONJ(H[0][2]);

      H[1][2].re = (H[1][2].re + H[2][1].re)/2.;
      H[1][2].im = (H[1][2].im - H[2][1].im)/2.;

      H[2][1] =  qcd_CONJ(H[1][2]);
      /*
        If H^2 is alread diagonal skip diagonalization and
        calculate U directly
      */
      sum=qcd_NORM(H[0][1])+qcd_NORM(H[0][2])+qcd_NORM(H[1][2]);

      if(sum<=1e-08)
      {
         e[0]=1./sqrt(H[0][0].re);
         e[1]=1./sqrt(H[1][1].re);
         e[2]=1./sqrt(H[2][2].re);

         U[0][0] = (qcd_complex_16) { M[0][0].re*e[0], M[0][0].im*e[0] };
         U[0][1] = (qcd_complex_16) { M[0][1].re*e[0], M[0][1].im*e[0] };
         U[0][2] = (qcd_complex_16) { M[0][2].re*e[0], M[0][2].im*e[0] };

         U[1][0] = (qcd_complex_16) { M[1][0].re*e[1], M[1][0].im*e[1] };
         U[1][1] = (qcd_complex_16) { M[1][1].re*e[1], M[1][1].im*e[1] };
         U[1][2] = (qcd_complex_16) { M[1][2].re*e[1], M[1][2].im*e[1] };

         U[2][0] = (qcd_complex_16) { M[2][0].re*e[2], M[2][0].im*e[2] };
         U[2][1] = (qcd_complex_16) { M[2][1].re*e[2], M[2][1].im*e[2] };
         U[2][2] = (qcd_complex_16) { M[2][2].re*e[2], M[2][2].im*e[2] };

         for(c1=0; c1<3; c1++)
         for(c2=0; c2<3; c2++)
            gf->D[iv][mu][c1][c2] = U[c1][c2]; 
      }
      else
      {
         /*
           Make traceless to eliminate second order term in eigenvalue equation,
           i.e. eig^3 + a eig + b = 0, when H is traceless.
         */
         trace=(H[0][0].re+H[1][1].re+H[2][2].re)/3.;

         H[0][0].re-=trace;
         H[1][1].re-=trace;
         H[2][2].re-=trace;


         /*
           Solve for eigenvalues:
           e^3 - e (H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2) - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 + H11*|H23|^2 - H13*H12^* *H23^*,
           e^3 + a*e + b = 0,

           a = -(H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2)
           b = - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 - H22*|H23|^2 - H33*|H23|^2 - H13*H12^* *H23^*

           D=(-9b + sqrt(12a^3 + 81b^2))^(1/3)

           e = D/(18^(1/3)) - ((2/3)^(1/3))/D
           e = (1 + I sqrt(3))a / (D 12^(1/3)) - (1 - I sqrt(3)) D / (2 18^(1/3))
           e = (1 - I sqrt(3))a / (D 12^(1/3)) - (1 + I sqrt(3)) D / (2 18^(1/3))
         */
         a = -(H[2][2].re*H[2][2].re - H[0][0].re*H[1][1].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + qcd_CMULR(H[0][2],qcd_CONJ(H[0][2])) + qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])));

         b.re  = - H[0][0].re*H[1][1].re*H[2][2].re + H[2][2].re*qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULR(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULR(H[0][2],qcd_CONJ(H[0][2]));
         b.im  =                                      H[2][2].re*qcd_CMULI(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULI(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULI(H[0][2],qcd_CONJ(H[0][2]));

         b.re +=   H[0][0].re*qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULR(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));
         b.im +=   H[0][0].re*qcd_CMULI(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULI(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));

         w.re=qcd_CPOWR(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);
         w.im=qcd_CPOWI(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);

         D=qcd_CPOW(((qcd_complex_16){-9.*b.re + w.re,-9.*b.im + w.im}),1./3.); 

         e[0] = D.re/(ThirdRoot_18) - qcd_CDEVR(((qcd_complex_16){a*ThirdRoot_2_3,0}),D);
         e[1] = a*qcd_CDEVR(ThirdRootOne[0],((qcd_complex_16){D.re*ThirdRoot_12,D.im*ThirdRoot_12})) - qcd_CMULR(ThirdRootOne[1],D)/(ThirdRoot_18*2.);
         e[2] = -e[0]-e[1];

         e[0]+= trace;
         e[1]+= trace;
         e[2]+= trace;

         H[0][0].re+=trace;
         H[1][1].re+=trace;
         H[2][2].re+=trace;

         /*
           Eigenvectors:
           v[0] = -(e H31 - H31 H22 + H21 H32) v[2] / Denom
           v[1] = -(H31 H12 - e H32 - H11 H32) v[2] / Denom
           v[2] =  (-e^2) + e H11 + |H12|^2 + e H22 - H11 H22
         */

         v[0][0].re = -(e[0]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[0][0].im = -(e[0]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[0][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[0]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[0][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[0]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[0][2].re =-e[0]*e[0] + e[0]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[0]*H[1][1].re - H[0][0].re*H[1][1].re;
         v[0][2].im = 0.;

         v[1][0].re = -(e[1]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
         v[1][0].im = -(e[1]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

         v[1][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[1]*H[2][1].re - H[0][0].re*H[2][1].re);
         v[1][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[1]*H[2][1].im - H[0][0].re*H[2][1].im);

         v[1][2].re =-e[1]*e[1] + e[1]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[1]*H[1][1].re - H[0][0].re*H[1][1].re;;
         v[1][2].im = 0.;

         /*
           Assure eigenvectors orthonormality:
           norm =  inner product v1.v1
           w    = (inner product v1.v2)/norm
           v2   = w*v1
         */

         norm  = qcd_CMULR(v[0][0],qcd_CONJ(v[0][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[0][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[0][2]));
         w.re  = qcd_CMULR(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[1][2]));
         w.im  = qcd_CMULI(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULI(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULI(v[0][2],qcd_CONJ(v[1][2]));
         w.re /= norm;
         w.im /= norm;

         v[1][0].re-= qcd_CMULR(w,v[0][0]);
         v[1][0].im-= qcd_CMULI(w,v[0][0]);

         v[1][1].re-= qcd_CMULR(w,v[0][1]);
         v[1][1].im-= qcd_CMULI(w,v[0][1]);

         v[1][2].re-= qcd_CMULR(w,v[0][2]);
         v[1][2].im-= qcd_CMULI(w,v[0][2]);

         norm=1./sqrt(norm);

         /*
           Normalize first and second eigenvector:
         */

         v[0][0].re*= norm;
         v[0][0].im*= norm;

         v[0][1].re*= norm;
         v[0][1].im*= norm;

         v[0][2].re*= norm;
         v[0][2].im*= norm;


         norm = qcd_CMULR(v[1][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[1][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[1][2],qcd_CONJ(v[1][2]));

         norm=1./sqrt(norm);

         v[1][0].re*= norm;
         v[1][0].im*= norm;

         v[1][1].re*= norm;
         v[1][1].im*= norm;

         v[1][2].re*= norm;
         v[1][2].im*= norm;

         /*
           v3 = v1 x v2
         */


         v[2][0].re =  qcd_CMULR(v[0][1],v[1][2]) - qcd_CMULR(v[0][2],v[1][1]);
         v[2][0].im = -qcd_CMULI(v[0][1],v[1][2]) + qcd_CMULI(v[0][2],v[1][1]);

         v[2][1].re = -qcd_CMULR(v[0][0],v[1][2]) + qcd_CMULR(v[0][2],v[1][0]);
         v[2][1].im = +qcd_CMULI(v[0][0],v[1][2]) - qcd_CMULI(v[0][2],v[1][0]);

         v[2][2].re =  qcd_CMULR(v[0][0],v[1][1]) - qcd_CMULR(v[0][1],v[1][0]);
         v[2][2].im = -qcd_CMULI(v[0][0],v[1][1]) + qcd_CMULI(v[0][1],v[1][0]);

         de     =               e[0]*e[1] +   e[1]*e[2] +   e[2]*e[0];
         /*
         cor[0] = tan(phase) * (e[0]*e[1] - 2*e[1]*e[2] +   e[2]*e[0])/de;
         cor[1] = tan(phase) * (e[0]*e[1] +   e[1]*e[2] - 2*e[2]*e[0])/de;
         cor[2] = - cor[0] - cor[1];
         */
         //to be compatible with Grenoble & Paris, don't apply corrections
         cor[0]=0;
         cor[1]=0;
         cor[2]=0;

         de = 1./sqrt(e[0]);
         b.re = de*cos(phase-cor[0]);
         b.im =-de*sin(phase-cor[0]);
         vr[0][0] = qcd_CMUL(b,v[0][0]);
         vr[0][1] = qcd_CMUL(b,v[0][1]);
         vr[0][2] = qcd_CMUL(b,v[0][2]);

         de = 1./sqrt(e[1]);
         b.re = de*cos(phase-cor[1]);
         b.im =-de*sin(phase-cor[1]);

         vr[1][0] = qcd_CMUL(b,v[1][0]);
         vr[1][1] = qcd_CMUL(b,v[1][1]);
         vr[1][2] = qcd_CMUL(b,v[1][2]);

         de = 1./sqrt(e[2]);
         b.re = de*cos(phase-cor[2]);
         b.im =-de*sin(phase-cor[2]);

         vr[2][0] = qcd_CMUL(b,v[2][0]);
         vr[2][1] = qcd_CMUL(b,v[2][1]);
         vr[2][2] = qcd_CMUL(b,v[2][2]);


         H[0][0].re= qcd_CMULR(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].re= qcd_CMULR(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].re= qcd_CMULR(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[2][2])) ;

         H[0][0].im= qcd_CMULI(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[0][2])) ;
         H[0][1].im= qcd_CMULI(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[1][2])) ;
         H[0][2].im= qcd_CMULI(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[2][2])) ;


         H[1][0].re= qcd_CMULR(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].re= qcd_CMULR(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].re= qcd_CMULR(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[2][2])) ;

         H[1][0].im= qcd_CMULI(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[0][2])) ;
         H[1][1].im= qcd_CMULI(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[1][2])) ;
         H[1][2].im= qcd_CMULI(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[2][2])) ;


         H[2][0].re= qcd_CMULR(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].re= qcd_CMULR(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].re= qcd_CMULR(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[2][2])) ;

         H[2][0].im= qcd_CMULI(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[0][2])) ;
         H[2][1].im= qcd_CMULI(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[1][2])) ;
         H[2][2].im= qcd_CMULI(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[2][2])) ;

         U[0][0].re= qcd_CMULR(H[0][0],vr[0][0]) +  qcd_CMULR(H[0][1],vr[1][0])  +  qcd_CMULR(H[0][2],vr[2][0]) ;
         U[0][1].re= qcd_CMULR(H[0][0],vr[0][1]) +  qcd_CMULR(H[0][1],vr[1][1])  +  qcd_CMULR(H[0][2],vr[2][1]) ;
         U[0][2].re= qcd_CMULR(H[0][0],vr[0][2]) +  qcd_CMULR(H[0][1],vr[1][2])  +  qcd_CMULR(H[0][2],vr[2][2]) ;

         U[0][0].im= qcd_CMULI(H[0][0],vr[0][0]) +  qcd_CMULI(H[0][1],vr[1][0])  +  qcd_CMULI(H[0][2],vr[2][0]) ;
         U[0][1].im= qcd_CMULI(H[0][0],vr[0][1]) +  qcd_CMULI(H[0][1],vr[1][1])  +  qcd_CMULI(H[0][2],vr[2][1]) ;
         U[0][2].im= qcd_CMULI(H[0][0],vr[0][2]) +  qcd_CMULI(H[0][1],vr[1][2])  +  qcd_CMULI(H[0][2],vr[2][2]) ;


         U[1][0].re= qcd_CMULR(H[1][0],vr[0][0]) +  qcd_CMULR(H[1][1],vr[1][0])  +  qcd_CMULR(H[1][2],vr[2][0]) ;
         U[1][1].re= qcd_CMULR(H[1][0],vr[0][1]) +  qcd_CMULR(H[1][1],vr[1][1])  +  qcd_CMULR(H[1][2],vr[2][1]) ;
         U[1][2].re= qcd_CMULR(H[1][0],vr[0][2]) +  qcd_CMULR(H[1][1],vr[1][2])  +  qcd_CMULR(H[1][2],vr[2][2]) ;

         U[1][0].im= qcd_CMULI(H[1][0],vr[0][0]) +  qcd_CMULI(H[1][1],vr[1][0])  +  qcd_CMULI(H[1][2],vr[2][0]) ;
         U[1][1].im= qcd_CMULI(H[1][0],vr[0][1]) +  qcd_CMULI(H[1][1],vr[1][1])  +  qcd_CMULI(H[1][2],vr[2][1]) ;
         U[1][2].im= qcd_CMULI(H[1][0],vr[0][2]) +  qcd_CMULI(H[1][1],vr[1][2])  +  qcd_CMULI(H[1][2],vr[2][2]) ;


         U[2][0].re= qcd_CMULR(H[2][0],vr[0][0]) +  qcd_CMULR(H[2][1],vr[1][0])  +  qcd_CMULR(H[2][2],vr[2][0]) ;
         U[2][1].re= qcd_CMULR(H[2][0],vr[0][1]) +  qcd_CMULR(H[2][1],vr[1][1])  +  qcd_CMULR(H[2][2],vr[2][1]) ;
         U[2][2].re= qcd_CMULR(H[2][0],vr[0][2]) +  qcd_CMULR(H[2][1],vr[1][2])  +  qcd_CMULR(H[2][2],vr[2][2]) ;

         U[2][0].im= qcd_CMULI(H[2][0],vr[0][0]) +  qcd_CMULI(H[2][1],vr[1][0])  +  qcd_CMULI(H[2][2],vr[2][0]) ;
         U[2][1].im= qcd_CMULI(H[2][0],vr[0][1]) +  qcd_CMULI(H[2][1],vr[1][1])  +  qcd_CMULI(H[2][2],vr[2][1]) ;
         U[2][2].im= qcd_CMULI(H[2][0],vr[0][2]) +  qcd_CMULI(H[2][1],vr[1][2])  +  qcd_CMULI(H[2][2],vr[2][2]) ;


         /*
           w    = inner product: col1.col2
           norm = inner product: col1.col1
         */

         norm  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][0])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][0])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][0]));
         w.re  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][1]));
         w.im  = qcd_CMULI(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULI(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULI(U[2][0],qcd_CONJ(U[2][1]));
         w.re /= norm;
         w.im /= norm;


         U[0][1].re-=qcd_CMULR(w,U[0][0]);
         U[0][1].im-=qcd_CMULI(w,U[0][0]);

         U[1][1].re-=qcd_CMULR(w,U[1][0]);
         U[1][1].im-=qcd_CMULI(w,U[1][0]);

         U[2][1].re-=qcd_CMULR(w,U[2][0]);
         U[2][1].im-=qcd_CMULI(w,U[2][0]);

         norm = 1./sqrt(norm);

         U[0][0].re*= norm;
         U[0][0].im*= norm;
         U[1][0].re*= norm;
         U[1][0].im*= norm;
         U[2][0].re*= norm;
         U[2][0].im*= norm;

         norm = qcd_CMULR(U[0][1],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][1],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][1],qcd_CONJ(U[2][1]));
         norm = 1./sqrt(norm);

         U[0][1].re*= norm;
         U[0][1].im*= norm;
         U[1][1].re*= norm;
         U[1][1].im*= norm;
         U[2][1].re*= norm;
         U[2][1].im*= norm;

         /*
           col3 = col1 x col2
         */
         U[0][2].re =  qcd_CMULR(U[1][0],U[2][1]) - qcd_CMULR(U[2][0],U[1][1]);
         U[0][2].im = -qcd_CMULI(U[1][0],U[2][1]) + qcd_CMULI(U[2][0],U[1][1]);

         U[1][2].re = -qcd_CMULR(U[0][0],U[2][1]) + qcd_CMULR(U[2][0],U[0][1]);
         U[1][2].im =  qcd_CMULI(U[0][0],U[2][1]) - qcd_CMULI(U[2][0],U[0][1]);

         U[2][2].re =  qcd_CMULR(U[0][0],U[1][1]) - qcd_CMULR(U[1][0],U[0][1]);
         U[2][2].im = -qcd_CMULI(U[0][0],U[1][1]) + qcd_CMULI(U[1][0],U[0][1]);

         for(c1=0; c1<3; c1++)
         for(c2=0; c2<3; c2++)
            gf->D[iv][mu][c1][c2] = U[c1][c2]; 
      }
   }//end volume-mu-loop
   return;
}//end qcd_projectSU33d


/*
  Projects the a single 3x3 matrix to SU(3). Written by Giannis Koutsou. Adapted to qcd_lib format by T.K.
  If M is an arbitary 3x3 matrix it can be factorised as the product of a hermitian and a unitary
  matrix, M = U H. We will calculate U as the projection of M to SU(3). Now M^\dagger M = H^2. By
  going to diagonal basis for H^2 we calculate H^-1 = 1./Sqrt(H^2) = v 1./Sqrt(L) v^\dagger where v
  is the 3x3 matrix whoes columns are the eigenvectors of H^2 and L is the diagonal representation
  of H^2. Thus: U = M H^-1
  [1] Phys. Lett. B307, 375-382, (1993)
*/
void qcd_projectSU33x3(qcd_complex_16 g[3][3])
{
   qcd_uint_4 iv,mu,c1,c2;
   qcd_complex_16 H[3][3],detM,U[3][3],M[3][3];
   qcd_complex_16 b,D,v[3][3],vr[3][3];
   qcd_real_8 a,ThirdRoot_18,ThirdRoot_12,ThirdRoot_2_3;
   qcd_real_8 trace,e[3],de,cor[3];
   qcd_complex_16 w,ThirdRootOne[2];
   qcd_real_8  sum;
   qcd_real_8  norm,phase;

   /*
     Constants Used:
   */
   ThirdRootOne[0].re= 1;
   ThirdRootOne[0].im= sqrt(3.);

   ThirdRootOne[1].re= 1;
   ThirdRootOne[1].im=-sqrt(3.);

   ThirdRoot_12 =pow(12.,1./3.);
   ThirdRoot_18 =pow(18.,1./3.);
   ThirdRoot_2_3=pow((2./3.),1./3.);

   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
      M[c1][c2] = g[c1][c2];


   detM = qcd_CADD(qcd_CADD( qcd_CMUL(M[0][0],qcd_CMUL(M[1][1],M[2][2])),
                             qcd_CMUL(M[0][1],qcd_CMUL(M[1][2],M[2][0]))),
                             qcd_CMUL(M[0][2],qcd_CMUL(M[1][0],M[2][1])));

   detM = qcd_CSUB(detM,
          qcd_CADD(qcd_CADD( qcd_CMUL(M[0][2],qcd_CMUL(M[1][1],M[2][0])),
                             qcd_CMUL(M[0][0],qcd_CMUL(M[1][2],M[2][1]))),
                             qcd_CMUL(M[0][1],qcd_CMUL(M[1][0],M[2][2]))));

   phase = qcd_ARG(detM)/3.;

   H[0][0].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][0]);
   H[0][1].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][1]);
   H[0][2].re= qcd_CMULR(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][0]),M[2][2]);
   H[0][0].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][0]);
   H[0][1].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][1]);
   H[0][2].im= qcd_CMULI(qcd_CONJ(M[0][0]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][0]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][0]),M[2][2]);

   H[1][0].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][0]);
   H[1][1].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][1]);
   H[1][2].re= qcd_CMULR(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][1]),M[2][2]);
   H[1][0].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][0]);
   H[1][1].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][1]);
   H[1][2].im= qcd_CMULI(qcd_CONJ(M[0][1]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][1]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][1]),M[2][2]);

   H[2][0].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][0]);
   H[2][1].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][1]);
   H[2][2].re= qcd_CMULR(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULR(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULR(qcd_CONJ(M[2][2]),M[2][2]);
   H[2][0].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][0]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][0]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][0]);
   H[2][1].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][1]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][1]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][1]);
   H[2][2].im= qcd_CMULI(qcd_CONJ(M[0][2]),M[0][2]) +  qcd_CMULI(qcd_CONJ(M[1][2]),M[1][2]) +  qcd_CMULI(qcd_CONJ(M[2][2]),M[2][2]);

   /*
     Assure Hermiticity:
   */
   H[0][1].re = (H[0][1].re + H[1][0].re)/2.;
   H[0][1].im = (H[0][1].im - H[1][0].im)/2.;

   H[1][0] =  qcd_CONJ(H[0][1]);

   H[0][2].re = (H[0][2].re + H[2][0].re)/2.;
   H[0][2].im = (H[0][2].im - H[2][0].im)/2.;

   H[2][0] =  qcd_CONJ(H[0][2]);

   H[1][2].re = (H[1][2].re + H[2][1].re)/2.;
   H[1][2].im = (H[1][2].im - H[2][1].im)/2.;

   H[2][1] =  qcd_CONJ(H[1][2]);
   /*
     If H^2 is alread diagonal skip diagonalization and
     calculate U directly
   */
   sum=qcd_NORM(H[0][1])+qcd_NORM(H[0][2])+qcd_NORM(H[1][2]);

   if(sum<=1e-12)
   {
      e[0]=1./sqrt(H[0][0].re);
      e[1]=1./sqrt(H[1][1].re);
      e[2]=1./sqrt(H[2][2].re);

      U[0][0] = (qcd_complex_16) { M[0][0].re*e[0], M[0][0].im*e[0] };
      U[0][1] = (qcd_complex_16) { M[0][1].re*e[0], M[0][1].im*e[0] };
      U[0][2] = (qcd_complex_16) { M[0][2].re*e[0], M[0][2].im*e[0] };

      U[1][0] = (qcd_complex_16) { M[1][0].re*e[1], M[1][0].im*e[1] };
      U[1][1] = (qcd_complex_16) { M[1][1].re*e[1], M[1][1].im*e[1] };
      U[1][2] = (qcd_complex_16) { M[1][2].re*e[1], M[1][2].im*e[1] };

      U[2][0] = (qcd_complex_16) { M[2][0].re*e[2], M[2][0].im*e[2] };
      U[2][1] = (qcd_complex_16) { M[2][1].re*e[2], M[2][1].im*e[2] };
      U[2][2] = (qcd_complex_16) { M[2][2].re*e[2], M[2][2].im*e[2] };

      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
         g[c1][c2] = U[c1][c2]; 
   }
   else
   {
      /*
        Make traceless to eliminate second order term in eigenvalue equation,
        i.e. eig^3 + a eig + b = 0, when H is traceless.
      */
      trace=(H[0][0].re+H[1][1].re+H[2][2].re)/3.;

      H[0][0].re-=trace;
      H[1][1].re-=trace;
      H[2][2].re-=trace;


      /*
        Solve for eigenvalues:
        e^3 - e (H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2) - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 + H11*|H23|^2 - H13*H12^* *H23^*,
        e^3 + a*e + b = 0,

        a = -(H33^2 - H11*H22 + |H12|^2 + |H13|^2 + |H23|^2)
        b = - H11*H22*H33 + H33*|H12|^2 - H12*H23*H13^* + H22*|H13|^2 - H22*|H23|^2 - H33*|H23|^2 - H13*H12^* *H23^*

        D=(-9b + sqrt(12a^3 + 81b^2))^(1/3)

        e = D/(18^(1/3)) - ((2/3)^(1/3))/D
        e = (1 + I sqrt(3))a / (D 12^(1/3)) - (1 - I sqrt(3)) D / (2 18^(1/3))
        e = (1 - I sqrt(3))a / (D 12^(1/3)) - (1 + I sqrt(3)) D / (2 18^(1/3))
      */
      a = -(H[2][2].re*H[2][2].re - H[0][0].re*H[1][1].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + qcd_CMULR(H[0][2],qcd_CONJ(H[0][2])) + qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])));

      b.re  = - H[0][0].re*H[1][1].re*H[2][2].re + H[2][2].re*qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULR(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULR(H[0][2],qcd_CONJ(H[0][2]));
      b.im  =                                      H[2][2].re*qcd_CMULI(H[0][1],qcd_CONJ(H[0][1])) - qcd_CMULI(H[0][1],qcd_CMUL(H[1][2],qcd_CONJ(H[0][2]))) + H[1][1].re*qcd_CMULI(H[0][2],qcd_CONJ(H[0][2]));

      b.re +=   H[0][0].re*qcd_CMULR(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULR(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));
      b.im +=   H[0][0].re*qcd_CMULI(H[1][2],qcd_CONJ(H[1][2])) - qcd_CMULI(H[0][2],qcd_CMUL(qcd_CONJ(H[0][1]),qcd_CONJ(H[1][2])));

      w.re=qcd_CPOWR(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);
      w.im=qcd_CPOWI(((qcd_complex_16){12.*a*a*a + 81.*qcd_CMULR(b,b), 81.*qcd_CMULI(b,b)}),0.5);

      D=qcd_CPOW(((qcd_complex_16){-9.*b.re + w.re,-9.*b.im + w.im}),1./3.); 

      e[0] = D.re/(ThirdRoot_18) - qcd_CDEVR(((qcd_complex_16){a*ThirdRoot_2_3,0}),D);
      e[1] = a*qcd_CDEVR(ThirdRootOne[0],((qcd_complex_16){D.re*ThirdRoot_12,D.im*ThirdRoot_12})) - qcd_CMULR(ThirdRootOne[1],D)/(ThirdRoot_18*2.);
      e[2] = -e[0]-e[1];

      e[0]+= trace;
      e[1]+= trace;
      e[2]+= trace;

      H[0][0].re+=trace;
      H[1][1].re+=trace;
      H[2][2].re+=trace;

      /*
        Eigenvectors:
        v[0] = -(e H31 - H31 H22 + H21 H32) v[2] / Denom
        v[1] = -(H31 H12 - e H32 - H11 H32) v[2] / Denom
        v[2] =  (-e^2) + e H11 + |H12|^2 + e H22 - H11 H22
      */

      v[0][0].re = -(e[0]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
      v[0][0].im = -(e[0]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

      v[0][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[0]*H[2][1].re - H[0][0].re*H[2][1].re);
      v[0][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[0]*H[2][1].im - H[0][0].re*H[2][1].im);

      v[0][2].re =-e[0]*e[0] + e[0]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[0]*H[1][1].re - H[0][0].re*H[1][1].re;
      v[0][2].im = 0.;

      v[1][0].re = -(e[1]*H[2][0].re - H[2][0].re*H[1][1].re + qcd_CMULR(H[1][0],H[2][1]));
      v[1][0].im = -(e[1]*H[2][0].im - H[2][0].im*H[1][1].re + qcd_CMULI(H[1][0],H[2][1]));

      v[1][1].re = -(qcd_CMULR(H[2][0],H[0][1]) + e[1]*H[2][1].re - H[0][0].re*H[2][1].re);
      v[1][1].im = -(qcd_CMULI(H[2][0],H[0][1]) + e[1]*H[2][1].im - H[0][0].re*H[2][1].im);

      v[1][2].re =-e[1]*e[1] + e[1]*H[0][0].re + qcd_CMULR(H[0][1],qcd_CONJ(H[0][1])) + e[1]*H[1][1].re - H[0][0].re*H[1][1].re;;
      v[1][2].im = 0.;

      /*
        Assure eigenvectors orthonormality:
        norm =  inner product v1.v1
        w    = (inner product v1.v2)/norm
        v2   = w*v1
      */

      norm  = qcd_CMULR(v[0][0],qcd_CONJ(v[0][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[0][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[0][2]));
      w.re  = qcd_CMULR(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[0][2],qcd_CONJ(v[1][2]));
      w.im  = qcd_CMULI(v[0][0],qcd_CONJ(v[1][0])) + qcd_CMULI(v[0][1],qcd_CONJ(v[1][1])) + qcd_CMULI(v[0][2],qcd_CONJ(v[1][2]));
      w.re /= norm;
      w.im /= norm;

      v[1][0].re-= qcd_CMULR(w,v[0][0]);
      v[1][0].im-= qcd_CMULI(w,v[0][0]);

      v[1][1].re-= qcd_CMULR(w,v[0][1]);
      v[1][1].im-= qcd_CMULI(w,v[0][1]);

      v[1][2].re-= qcd_CMULR(w,v[0][2]);
      v[1][2].im-= qcd_CMULI(w,v[0][2]);

      norm=1./sqrt(norm);

      /*
        Normalize first and second eigenvector:
      */

      v[0][0].re*= norm;
      v[0][0].im*= norm;

      v[0][1].re*= norm;
      v[0][1].im*= norm;

      v[0][2].re*= norm;
      v[0][2].im*= norm;


      norm = qcd_CMULR(v[1][0],qcd_CONJ(v[1][0])) + qcd_CMULR(v[1][1],qcd_CONJ(v[1][1])) + qcd_CMULR(v[1][2],qcd_CONJ(v[1][2]));

      norm=1./sqrt(norm);

      v[1][0].re*= norm;
      v[1][0].im*= norm;

      v[1][1].re*= norm;
      v[1][1].im*= norm;

      v[1][2].re*= norm;
      v[1][2].im*= norm;

      /*
        v3 = v1 x v2
      */


      v[2][0].re =  qcd_CMULR(v[0][1],v[1][2]) - qcd_CMULR(v[0][2],v[1][1]);
      v[2][0].im = -qcd_CMULI(v[0][1],v[1][2]) + qcd_CMULI(v[0][2],v[1][1]);

      v[2][1].re = -qcd_CMULR(v[0][0],v[1][2]) + qcd_CMULR(v[0][2],v[1][0]);
      v[2][1].im = +qcd_CMULI(v[0][0],v[1][2]) - qcd_CMULI(v[0][2],v[1][0]);

      v[2][2].re =  qcd_CMULR(v[0][0],v[1][1]) - qcd_CMULR(v[0][1],v[1][0]);
      v[2][2].im = -qcd_CMULI(v[0][0],v[1][1]) + qcd_CMULI(v[0][1],v[1][0]);

      de     =               e[0]*e[1] +   e[1]*e[2] +   e[2]*e[0];
      /*
      cor[0] = tan(phase) * (e[0]*e[1] - 2*e[1]*e[2] +   e[2]*e[0])/de;
      cor[1] = tan(phase) * (e[0]*e[1] +   e[1]*e[2] - 2*e[2]*e[0])/de;
      cor[2] = - cor[0] - cor[1];
      */
      //to be compatible with Grenoble & Paris, don't apply corrections
      cor[0]=0;
      cor[1]=0;
      cor[2]=0;

      de = 1./sqrt(e[0]);
      b.re = de*cos(phase-cor[0]);
      b.im =-de*sin(phase-cor[0]);
      vr[0][0] = qcd_CMUL(b,v[0][0]);
      vr[0][1] = qcd_CMUL(b,v[0][1]);
      vr[0][2] = qcd_CMUL(b,v[0][2]);

      de = 1./sqrt(e[1]);
      b.re = de*cos(phase-cor[1]);
      b.im =-de*sin(phase-cor[1]);

      vr[1][0] = qcd_CMUL(b,v[1][0]);
      vr[1][1] = qcd_CMUL(b,v[1][1]);
      vr[1][2] = qcd_CMUL(b,v[1][2]);

      de = 1./sqrt(e[2]);
      b.re = de*cos(phase-cor[2]);
      b.im =-de*sin(phase-cor[2]);

      vr[2][0] = qcd_CMUL(b,v[2][0]);
      vr[2][1] = qcd_CMUL(b,v[2][1]);
      vr[2][2] = qcd_CMUL(b,v[2][2]);


      H[0][0].re= qcd_CMULR(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[0][2])) ;
      H[0][1].re= qcd_CMULR(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[1][2])) ;
      H[0][2].re= qcd_CMULR(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[0][2],qcd_CONJ(v[2][2])) ;

      H[0][0].im= qcd_CMULI(M[0][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[0][2])) ;
      H[0][1].im= qcd_CMULI(M[0][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[1][2])) ;
      H[0][2].im= qcd_CMULI(M[0][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[0][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[0][2],qcd_CONJ(v[2][2])) ;


      H[1][0].re= qcd_CMULR(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[0][2])) ;
      H[1][1].re= qcd_CMULR(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[1][2])) ;
      H[1][2].re= qcd_CMULR(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[1][2],qcd_CONJ(v[2][2])) ;

      H[1][0].im= qcd_CMULI(M[1][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[0][2])) ;
      H[1][1].im= qcd_CMULI(M[1][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[1][2])) ;
      H[1][2].im= qcd_CMULI(M[1][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[1][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[1][2],qcd_CONJ(v[2][2])) ;


      H[2][0].re= qcd_CMULR(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[0][2])) ;
      H[2][1].re= qcd_CMULR(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[1][2])) ;
      H[2][2].re= qcd_CMULR(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULR(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULR(M[2][2],qcd_CONJ(v[2][2])) ;

      H[2][0].im= qcd_CMULI(M[2][0],qcd_CONJ(v[0][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[0][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[0][2])) ;
      H[2][1].im= qcd_CMULI(M[2][0],qcd_CONJ(v[1][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[1][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[1][2])) ;
      H[2][2].im= qcd_CMULI(M[2][0],qcd_CONJ(v[2][0])) +  qcd_CMULI(M[2][1],qcd_CONJ(v[2][1]))  +  qcd_CMULI(M[2][2],qcd_CONJ(v[2][2])) ;

      U[0][0].re= qcd_CMULR(H[0][0],vr[0][0]) +  qcd_CMULR(H[0][1],vr[1][0])  +  qcd_CMULR(H[0][2],vr[2][0]) ;
      U[0][1].re= qcd_CMULR(H[0][0],vr[0][1]) +  qcd_CMULR(H[0][1],vr[1][1])  +  qcd_CMULR(H[0][2],vr[2][1]) ;
      U[0][2].re= qcd_CMULR(H[0][0],vr[0][2]) +  qcd_CMULR(H[0][1],vr[1][2])  +  qcd_CMULR(H[0][2],vr[2][2]) ;

      U[0][0].im= qcd_CMULI(H[0][0],vr[0][0]) +  qcd_CMULI(H[0][1],vr[1][0])  +  qcd_CMULI(H[0][2],vr[2][0]) ;
      U[0][1].im= qcd_CMULI(H[0][0],vr[0][1]) +  qcd_CMULI(H[0][1],vr[1][1])  +  qcd_CMULI(H[0][2],vr[2][1]) ;
      U[0][2].im= qcd_CMULI(H[0][0],vr[0][2]) +  qcd_CMULI(H[0][1],vr[1][2])  +  qcd_CMULI(H[0][2],vr[2][2]) ;


      U[1][0].re= qcd_CMULR(H[1][0],vr[0][0]) +  qcd_CMULR(H[1][1],vr[1][0])  +  qcd_CMULR(H[1][2],vr[2][0]) ;
      U[1][1].re= qcd_CMULR(H[1][0],vr[0][1]) +  qcd_CMULR(H[1][1],vr[1][1])  +  qcd_CMULR(H[1][2],vr[2][1]) ;
      U[1][2].re= qcd_CMULR(H[1][0],vr[0][2]) +  qcd_CMULR(H[1][1],vr[1][2])  +  qcd_CMULR(H[1][2],vr[2][2]) ;

      U[1][0].im= qcd_CMULI(H[1][0],vr[0][0]) +  qcd_CMULI(H[1][1],vr[1][0])  +  qcd_CMULI(H[1][2],vr[2][0]) ;
      U[1][1].im= qcd_CMULI(H[1][0],vr[0][1]) +  qcd_CMULI(H[1][1],vr[1][1])  +  qcd_CMULI(H[1][2],vr[2][1]) ;
      U[1][2].im= qcd_CMULI(H[1][0],vr[0][2]) +  qcd_CMULI(H[1][1],vr[1][2])  +  qcd_CMULI(H[1][2],vr[2][2]) ;


      U[2][0].re= qcd_CMULR(H[2][0],vr[0][0]) +  qcd_CMULR(H[2][1],vr[1][0])  +  qcd_CMULR(H[2][2],vr[2][0]) ;
      U[2][1].re= qcd_CMULR(H[2][0],vr[0][1]) +  qcd_CMULR(H[2][1],vr[1][1])  +  qcd_CMULR(H[2][2],vr[2][1]) ;
      U[2][2].re= qcd_CMULR(H[2][0],vr[0][2]) +  qcd_CMULR(H[2][1],vr[1][2])  +  qcd_CMULR(H[2][2],vr[2][2]) ;

      U[2][0].im= qcd_CMULI(H[2][0],vr[0][0]) +  qcd_CMULI(H[2][1],vr[1][0])  +  qcd_CMULI(H[2][2],vr[2][0]) ;
      U[2][1].im= qcd_CMULI(H[2][0],vr[0][1]) +  qcd_CMULI(H[2][1],vr[1][1])  +  qcd_CMULI(H[2][2],vr[2][1]) ;
      U[2][2].im= qcd_CMULI(H[2][0],vr[0][2]) +  qcd_CMULI(H[2][1],vr[1][2])  +  qcd_CMULI(H[2][2],vr[2][2]) ;


      /*
        w    = inner product: col1.col2
        norm = inner product: col1.col1
      */

      norm  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][0])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][0])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][0]));
      w.re  = qcd_CMULR(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][0],qcd_CONJ(U[2][1]));
      w.im  = qcd_CMULI(U[0][0],qcd_CONJ(U[0][1])) + qcd_CMULI(U[1][0],qcd_CONJ(U[1][1])) + qcd_CMULI(U[2][0],qcd_CONJ(U[2][1]));
      w.re /= norm;
      w.im /= norm;


      U[0][1].re-=qcd_CMULR(w,U[0][0]);
      U[0][1].im-=qcd_CMULI(w,U[0][0]);

      U[1][1].re-=qcd_CMULR(w,U[1][0]);
      U[1][1].im-=qcd_CMULI(w,U[1][0]);

      U[2][1].re-=qcd_CMULR(w,U[2][0]);
      U[2][1].im-=qcd_CMULI(w,U[2][0]);

      norm = 1./sqrt(norm);

      U[0][0].re*= norm;
      U[0][0].im*= norm;
      U[1][0].re*= norm;
      U[1][0].im*= norm;
      U[2][0].re*= norm;
      U[2][0].im*= norm;

      norm = qcd_CMULR(U[0][1],qcd_CONJ(U[0][1])) + qcd_CMULR(U[1][1],qcd_CONJ(U[1][1])) + qcd_CMULR(U[2][1],qcd_CONJ(U[2][1]));
      norm = 1./sqrt(norm);

      U[0][1].re*= norm;
      U[0][1].im*= norm;
      U[1][1].re*= norm;
      U[1][1].im*= norm;
      U[2][1].re*= norm;
      U[2][1].im*= norm;

      /*
        col3 = col1 x col2
      */
      U[0][2].re =  qcd_CMULR(U[1][0],U[2][1]) - qcd_CMULR(U[2][0],U[1][1]);
      U[0][2].im = -qcd_CMULI(U[1][0],U[2][1]) + qcd_CMULI(U[2][0],U[1][1]);

      U[1][2].re = -qcd_CMULR(U[0][0],U[2][1]) + qcd_CMULR(U[2][0],U[0][1]);
      U[1][2].im =  qcd_CMULI(U[0][0],U[2][1]) - qcd_CMULI(U[2][0],U[0][1]);

      U[2][2].re =  qcd_CMULR(U[0][0],U[1][1]) - qcd_CMULR(U[1][0],U[0][1]);
      U[2][2].im = -qcd_CMULI(U[0][0],U[1][1]) + qcd_CMULI(U[1][0],U[0][1]);

      for(c1=0; c1<3; c1++)
      for(c2=0; c2<3; c2++)
         g[c1][c2] = U[c1][c2];
   }
   return;
}//end qcd_projectSU33x3

void qcd_randomGaugeField(qcd_gaugeField *u)
{
   qcd_uint_4 x;
   qcd_uint_2 mu,c1,c2;
   for(x=0;x<u->geo->lV;x++)
   for(mu=0; mu<4; mu++)
   for(c1=0; c1<3; c1++)
   for(c2=0; c2<3; c2++)
   {
      u->D[x][mu][c1][c2] = (qcd_complex_16){drand48(),drand48()};
   }
   qcd_projectSU3(u);
}//end qcd_random
