/* qcd.h
 *
 * main header file for qcd-library
 * version for i386 machines
 *
 * Tomasz Korzec 2008
 ******************************************/

#ifndef _QCD_H

 
// some constants:
#define qcd_GF_LIME 0 

#define qcd_PROP_CMI 0
 
#define qcd_MAX_STRING_LENGTH 256 
 
// the following definitions are machine dependent: 
typedef float  qcd_real_4;              // 4 byte, single precision
typedef double qcd_real_8;              // 8 byte, double precision
typedef char   qcd_int_1;               // signed 1 byte integer
typedef short  qcd_int_2;               // signed 2 byte integer
typedef int    qcd_int_4;               // signed 4 byte integer
typedef long   qcd_int_8;               // signed 8 byte integer
typedef unsigned char   qcd_uint_1;     // unsigned 1 byte integer
typedef unsigned short  qcd_uint_2;     // unsigned 2 byte integer
typedef unsigned int    qcd_uint_4;     // unsigned 4 byte integer
typedef unsigned long   qcd_uint_8;     // unsigned 8 byte integer
//------------------------------------------------------------------




typedef struct
{
  qcd_real_8 re;
  qcd_real_8 im;
} qcd_complex_16; 


typedef struct 
{
  qcd_real_4 re;
  qcd_real_4 im;
} qcd_complex_8; 

typedef struct
{
   qcd_uint_8 V;                   // 4d volume
   qcd_uint_8 V3;                  // spatial volume
   qcd_uint_2 L[4];                // lattice dimensions t,x,y,z = 0,1,2,3
   qcd_uint_8 lV;                  // local volume
   qcd_uint_8 lV3;                 // local spatial volume
   qcd_uint_2 lL[4];               // dimensions of local lattice
   qcd_uint_2 Pos[4];              // position of local lattice
   qcd_real_8 theta[4];            // boundary conditions of fermion fields
                                   // psi(x+L[mu]) = exp(i*theta[mu])*psi(x)
   qcd_uint_4 myid;                // lexicographical index of position   
   qcd_uint_4 nproc;               // number of processes
   qcd_uint_4 Pplus[4];            // process IDs of neighbor blocks in + directions
   qcd_uint_4 Pminus[4];           // process IDs of neighbor blocks in - directions
   qcd_uint_8 (*plus)[4];          // neighborhood relations on local lattice
   qcd_uint_8 (*minus)[4];         //
   qcd_uint_8 startBminus[4];      // index where - boundaries start
   qcd_uint_8 startBplus[4];       // index where + boundaries start
   
   MPI_Request requests[32];       // needed in communication routines
   int numOfRequests;              // non 0 during ongoing communication
   MPI_Datatype stypeV[4];         // data-types that allow sends of strided data
   MPI_Datatype rtypeV[4];         //
   MPI_Datatype stypeU[4];         //
   MPI_Datatype rtypeU[4];         //
   
   qcd_uint_2 initialized;         // set to 1 during initialization 
} qcd_geometry;


typedef struct 
{
   qcd_complex_16 (*D)[4][3];         // data
   qcd_complex_16 (*Bplus[4])[4][3];  // boundary data in + direction
   qcd_complex_16 (*Bminus[4])[4][3]; // boundary data in - direction
   qcd_geometry *geo;                 // 
   qcd_uint_2 initialized;            // set to 1 during initialization
} qcd_vector;

typedef struct 
{
   qcd_complex_16 (*D)[4][4][3][3];         // data
   qcd_complex_16 (*Bplus[4])[4][4][3][3];  // boundary data in + direction
   qcd_complex_16 (*Bminus[4])[4][4][3][3]; // boundary data in - direction
   qcd_geometry *geo;                       // 
   qcd_uint_2 initialized;                  // set to 1 during initialization
} qcd_propagator;

typedef struct
{
   qcd_complex_16 (*D)[4][3][3];        // data
   qcd_complex_16 (*Bplus[4])[4][3][3]; // boundary data in + direction
   qcd_complex_16 (*Bminus[4])[4][3][3];// boundary data in - direction
   qcd_geometry *geo;
   qcd_uint_2 initialized;              // set to 1 during initialization
} qcd_gaugeField;



//////////////////////////////////////////////////////////////////////////////
// some pre-processor macros. Use carefully!!
#define qcd_LEXIC(t,x,y,z,L) ( (qcd_uint_8) ((t)+L[0]*((x)+L[1]*((y)+L[2]*(z)))) )
#define qcd_LEXIC0(x,y,z,L) ( (qcd_uint_8) ((x)+L[1]*((y)+L[2]*(z))) )
#define qcd_LEXIC1(t,y,z,L) ( (qcd_uint_8) ((t)+L[0]*((y)+L[2]*(z))) )
#define qcd_LEXIC2(t,x,z,L) ( (qcd_uint_8) ((t)+L[0]*((x)+L[1]*(z))) )
#define qcd_LEXIC3(t,x,y,L) ( (qcd_uint_8) ((t)+L[0]*((x)+L[1]*(y))) )

#define qcd_CONJ(x)    ( (qcd_complex_16) {x.re,-x.im} )
#define qcd_CMUL(x,y)  ( (qcd_complex_16) {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re } )
#define qcd_CMULR(x,y) ( (qcd_real_8) (x.re * y.re - x.im * y.im) )
#define qcd_CMULI(x,y) ( (qcd_real_8) (x.re * y.im + x.im * y.re) )

#define qcd_CADD(x,y)  ( (qcd_complex_16) {x.re+y.re, x.im+y.im})
#define qcd_CADDR(x,y) ( (qcd_real_8) (x.re+y.re))
#define qcd_CADDI(x,y) ( (qcd_real_8) (x.im+y.im))

#define qcd_CSUB(x,y)  ( (qcd_complex_16) {x.re-y.re, x.im-y.im})
#define qcd_CSUBR(x,y) ( (qcd_real_8) (x.re-y.re))
#define qcd_CSUBI(x,y) ( (qcd_real_8) (x.im-y.im))

#define qcd_CSCALE(x,a)  ( (qcd_complex_16) {x.re*(a), x.im*(a)})
#define qcd_CSCALER(x,a) ( (qcd_real_8) x.re*(a))
#define qcd_CSCALEI(x,a) ( (qcd_real_8) x.im*(a))

#define qcd_CDEV(x,y)  ( (qcd_complex_16) {(x.re * y.re + x.im * y.im)/(y.re*y.re + y.im*y.im), (x.im * y.re - x.re * y.im)/(y.re*y.re + y.im*y.im) } )
#define qcd_CDEVR(x,y) ( (qcd_real_8)     ((x.re * y.re + x.im * y.im)/(y.re*y.re + y.im*y.im) ))
#define qcd_CDEVI(x,y) ( (qcd_real_8)     ((x.im * y.re - x.re * y.im)/(y.re*y.re + y.im*y.im) ))

#define qcd_NORM(x)    ( (qcd_real_8) sqrt(x.re * x.re + x.im * x.im))
#define qcd_ARG(x)     ( (qcd_real_8) atan2(x.im,x.re))
#define qcd_CPOW(x,a)  ( (qcd_complex_16) {pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)), pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a))})
#define qcd_CPOWR(x,a) ( (qcd_real_8) pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)))
#define qcd_CPOWI(x,a) ( (qcd_real_8) pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a)))

#define _QCD_H 1
#endif
