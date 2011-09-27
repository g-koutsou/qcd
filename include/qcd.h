/* qcd.h
 *
 * main header file for qcd-library
 * version for i386 machines
 *
 * Tomasz Korzec 2008
 ******************************************/

#ifndef H_QCD
#define H_QCD 1

 
// some constants:

 
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
   qcd_uint_8 (*plus3d)[4];        // neighborhood relations on local lattice (timeslice)
   qcd_uint_8 (*minus3d)[4];       //
   qcd_uint_8 (*eo)[2];            // partitioning into even and odd 
   qcd_uint_8 startBminus[4];      // index where - boundaries start
   qcd_uint_8 startBplus[4];       // index where + boundaries start
   qcd_uint_8 startBminus3d[4];    // index where - boundaries start (timeslice)
   qcd_uint_8 startBplus3d[4];     // index where + boundaries start (timeslice)  
   qcd_uint_8 boundarySize;        // number of points in outer boundaries
   qcd_uint_8 boundarySize3d;      // number of points in outer boundaries (timeslice)
   qcd_uint_8 *edge;               // lexicographical indices of inner boundaries
   qcd_uint_8 edgePoints;          // number of inner boundary points
   qcd_uint_8 *edge0;              // lexicographical indices of inner b. of a time-slice
   qcd_uint_8 edge0Points;         // number of inner time-slice boundary points
   
   MPI_Request requests[32];       // needed in communication routines
   int numOfRequests;              // non 0 during ongoing communication
   MPI_Datatype stypeV[4];         // data-types that allow sends of strided data
   MPI_Datatype rtypeV[4];         //      vector-fields
   MPI_Datatype stypeU[4];         //      
   MPI_Datatype rtypeU[4];         //      gauge-fields
   MPI_Datatype stypeT[4];         //      
   MPI_Datatype rtypeT[4];         //      gauge-transformation-fields
   MPI_Datatype stypeP[4];         //
   MPI_Datatype rtypeP[4];         //      propagator-fields
   
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
   qcd_complex_16 (*D)[4][3];           // data
   qcd_complex_16 (*Bplus[4])[4][3];    // boundary data in + direction (only 1..3 are used)
   qcd_complex_16 (*Bminus[4])[4][3];   // boundary data in - direction
   qcd_geometry *geo;                   // 
   qcd_uint_2 initialized;              // set to 1 during initialization
} qcd_vector3d;


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

typedef struct
{
   qcd_complex_16 (*D)[3][3];        // data
   qcd_complex_16 (*Bplus[4])[3][3]; // boundary data in + direction
   qcd_complex_16 (*Bminus[4])[3][3];// boundary data in - direction
   qcd_geometry *geo;
   qcd_uint_2 initialized;              // set to 1 during initialization
} qcd_gaugeTransformation;


//////////////////////////////////////////////////////////////////////////////
// some pre-processor macros. Use carefully!!
#define qcd_LEXIC(t,x,y,z,L) ( (qcd_uint_8) ((t)+L[0]*((x)+L[1]*((y)+L[2]*(z)))) )
#define qcd_LEXIC0(x,y,z,L) ( (qcd_uint_8) ((x)+L[1]*((y)+L[2]*(z))) )
#define qcd_LEXIC1(t,y,z,L) ( (qcd_uint_8) ((t)+L[0]*((y)+L[2]*(z))) )
#define qcd_LEXIC2(t,x,z,L) ( (qcd_uint_8) ((t)+L[0]*((x)+L[1]*(z))) )
#define qcd_LEXIC3(t,x,y,L) ( (qcd_uint_8) ((t)+L[0]*((x)+L[1]*(y))) )
#define qcd_LEXIC01(y,z,L)  ( (qcd_uint_8) ((y)+L[2]*(z)) )
#define qcd_LEXIC02(x,z,L)  ( (qcd_uint_8) ((x)+L[1]*(z)) )
#define qcd_LEXIC03(x,y,L)  ( (qcd_uint_8) ((x)+L[1]*(y)) )



#include <qcd_gamma.h>
#include <qcd_communication.h>
#include <qcd_blas.h>
#include <qcd_io.h>
#include <qcd_smearing.h>
#include <qcd_gaugeFixing.h>
#include <qcd_wilson.h>
#include <qcd_observables.h>

/* prototypes for qcd_init.c*/
void qcd_antilexic(qcd_uint_2 x[], qcd_uint_8 l, const qcd_uint_2 dim[]);
int qcd_initGeometry(qcd_geometry *geo, const qcd_uint_2 L[], const qcd_uint_2 P[], const qcd_real_8 theta[], const int myid, const int numprocs);
void qcd_destroyGeometry(const qcd_geometry *geo);
int qcd_initVector(qcd_vector *vec, qcd_geometry *geo);
void qcd_destroyVector(qcd_vector *vec);
int qcd_initVector3d(qcd_vector3d *vec, qcd_geometry *geo);
void qcd_destroyVector3d(qcd_vector3d *vec);
int qcd_initPropagator(qcd_propagator *prp, qcd_geometry *geo);
void qcd_destroyPropagator(qcd_propagator *prp);
int qcd_initGaugeField(qcd_gaugeField *u, qcd_geometry *geo);
void qcd_destroyGaugeField(qcd_gaugeField *u);
int qcd_initGaugeTransformation(qcd_gaugeTransformation *u, qcd_geometry *geo);
void qcd_destroyGaugeTransformation(qcd_gaugeTransformation *u);

#endif
