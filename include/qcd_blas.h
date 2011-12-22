/* qcd_blas.h
 *
 * header file for qcd_blas.c
 *
 * Tomasz Korzec 2009
 *******************************************/
 
#ifndef H_QCD_BLAS
#define H_QCD_BLAS 1 
 
/* preprocessor macros, use carefully */
#define qcd_CONJ(x)    ( (qcd_complex_16) {x.re,-x.im} )
#define qcd_CMUL(x,y)  ( (qcd_complex_16) {x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re } )
#define qcd_CMULR(x,y) ( (qcd_real_8) (x.re * y.re - x.im * y.im) )
#define qcd_CMULI(x,y) ( (qcd_real_8) (x.re * y.im + x.im * y.re) )
#define qcd_CADJOINTMUL(x,y)  ( (qcd_complex_16) {x.re * y.re + x.im * y.im, x.re * y.im - x.im * y.re } )

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
#define qcd_NORMSQUARED(x) ( (qcd_real_8) (x.re * x.re + x.im * x.im))
#define qcd_ARG(x)     ( (qcd_real_8) atan2(x.im,x.re))
#define qcd_CPOW(x,a)  ( (qcd_complex_16) {pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)), pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a))})
#define qcd_CPOWR(x,a) ( (qcd_real_8) pow(qcd_NORM(x),(a))*cos(qcd_ARG(x)*(a)))
#define qcd_CPOWI(x,a) ( (qcd_real_8) pow(qcd_NORM(x),(a))*sin(qcd_ARG(x)*(a))) 
 
 
/* a = b*c, where a,b,c are 3x3 matrices */ 
#define qcd_MUL3x3(a,b,c) ({ a[0][0].re =  b[0][0].re*c[0][0].re - b[0][0].im*c[0][0].im \
                                          +b[0][1].re*c[1][0].re - b[0][1].im*c[1][0].im \
                                          +b[0][2].re*c[2][0].re - b[0][2].im*c[2][0].im;\
                             a[0][0].im =  b[0][0].re*c[0][0].im + b[0][0].im*c[0][0].re \
                                          +b[0][1].re*c[1][0].im + b[0][1].im*c[1][0].re \
                                          +b[0][2].re*c[2][0].im + b[0][2].im*c[2][0].re;\
                             a[0][1].re =  b[0][0].re*c[0][1].re - b[0][0].im*c[0][1].im \
                                          +b[0][1].re*c[1][1].re - b[0][1].im*c[1][1].im \
                                          +b[0][2].re*c[2][1].re - b[0][2].im*c[2][1].im;\
                             a[0][1].im =  b[0][0].re*c[0][1].im + b[0][0].im*c[0][1].re \
                                          +b[0][1].re*c[1][1].im + b[0][1].im*c[1][1].re \
                                          +b[0][2].re*c[2][1].im + b[0][2].im*c[2][1].re;\
                             a[0][2].re =  b[0][0].re*c[0][2].re - b[0][0].im*c[0][2].im \
                                          +b[0][1].re*c[1][2].re - b[0][1].im*c[1][2].im \
                                          +b[0][2].re*c[2][2].re - b[0][2].im*c[2][2].im;\
                             a[0][2].im =  b[0][0].re*c[0][2].im + b[0][0].im*c[0][2].re \
                                          +b[0][1].re*c[1][2].im + b[0][1].im*c[1][2].re \
                                          +b[0][2].re*c[2][2].im + b[0][2].im*c[2][2].re;\
                                                                                         \
                             a[1][0].re =  b[1][0].re*c[0][0].re - b[1][0].im*c[0][0].im \
                                          +b[1][1].re*c[1][0].re - b[1][1].im*c[1][0].im \
                                          +b[1][2].re*c[2][0].re - b[1][2].im*c[2][0].im;\
                             a[1][0].im =  b[1][0].re*c[0][0].im + b[1][0].im*c[0][0].re \
                                          +b[1][1].re*c[1][0].im + b[1][1].im*c[1][0].re \
                                          +b[1][2].re*c[2][0].im + b[1][2].im*c[2][0].re;\
                             a[1][1].re =  b[1][0].re*c[0][1].re - b[1][0].im*c[0][1].im \
                                          +b[1][1].re*c[1][1].re - b[1][1].im*c[1][1].im \
                                          +b[1][2].re*c[2][1].re - b[1][2].im*c[2][1].im;\
                             a[1][1].im =  b[1][0].re*c[0][1].im + b[1][0].im*c[0][1].re \
                                          +b[1][1].re*c[1][1].im + b[1][1].im*c[1][1].re \
                                          +b[1][2].re*c[2][1].im + b[1][2].im*c[2][1].re;\
                             a[1][2].re =  b[1][0].re*c[0][2].re - b[1][0].im*c[0][2].im \
                                          +b[1][1].re*c[1][2].re - b[1][1].im*c[1][2].im \
                                          +b[1][2].re*c[2][2].re - b[1][2].im*c[2][2].im;\
                             a[1][2].im =  b[1][0].re*c[0][2].im + b[1][0].im*c[0][2].re \
                                          +b[1][1].re*c[1][2].im + b[1][1].im*c[1][2].re \
                                          +b[1][2].re*c[2][2].im + b[1][2].im*c[2][2].re;\
                                                                                         \
                             a[2][0].re =  b[2][0].re*c[0][0].re - b[2][0].im*c[0][0].im \
                                          +b[2][1].re*c[1][0].re - b[2][1].im*c[1][0].im \
                                          +b[2][2].re*c[2][0].re - b[2][2].im*c[2][0].im;\
                             a[2][0].im =  b[2][0].re*c[0][0].im + b[2][0].im*c[0][0].re \
                                          +b[2][1].re*c[1][0].im + b[2][1].im*c[1][0].re \
                                          +b[2][2].re*c[2][0].im + b[2][2].im*c[2][0].re;\
                             a[2][1].re =  b[2][0].re*c[0][1].re - b[2][0].im*c[0][1].im \
                                          +b[2][1].re*c[1][1].re - b[2][1].im*c[1][1].im \
                                          +b[2][2].re*c[2][1].re - b[2][2].im*c[2][1].im;\
                             a[2][1].im =  b[2][0].re*c[0][1].im + b[2][0].im*c[0][1].re \
                                          +b[2][1].re*c[1][1].im + b[2][1].im*c[1][1].re \
                                          +b[2][2].re*c[2][1].im + b[2][2].im*c[2][1].re;\
                             a[2][2].re =  b[2][0].re*c[0][2].re - b[2][0].im*c[0][2].im \
                                          +b[2][1].re*c[1][2].re - b[2][1].im*c[1][2].im \
                                          +b[2][2].re*c[2][2].re - b[2][2].im*c[2][2].im;\
                             a[2][2].im =  b[2][0].re*c[0][2].im + b[2][0].im*c[0][2].re \
                                          +b[2][1].re*c[1][2].im + b[2][1].im*c[1][2].re \
                                          +b[2][2].re*c[2][2].im + b[2][2].im*c[2][2].re;\
})
 

/* a = b* adjoint(c), where a,b,c are 3x3 matrices */
#define qcd_MULADJOINT3x3(a,b,c) ({\
                             a[0][0].re =  b[0][0].re*c[0][0].re + b[0][0].im*c[0][0].im \
                                          +b[0][1].re*c[0][1].re + b[0][1].im*c[0][1].im \
                                          +b[0][2].re*c[0][2].re + b[0][2].im*c[0][2].im;\
                             a[0][0].im = -b[0][0].re*c[0][0].im + b[0][0].im*c[0][0].re \
                                          -b[0][1].re*c[0][1].im + b[0][1].im*c[0][1].re \
                                          -b[0][2].re*c[0][2].im + b[0][2].im*c[0][2].re;\
                             a[0][1].re =  b[0][0].re*c[1][0].re + b[0][0].im*c[1][0].im \
                                          +b[0][1].re*c[1][1].re + b[0][1].im*c[1][1].im \
                                          +b[0][2].re*c[1][2].re + b[0][2].im*c[1][2].im;\
                             a[0][1].im = -b[0][0].re*c[1][0].im + b[0][0].im*c[1][0].re \
                                          -b[0][1].re*c[1][1].im + b[0][1].im*c[1][1].re \
                                          -b[0][2].re*c[1][2].im + b[0][2].im*c[1][2].re;\
                             a[0][2].re =  b[0][0].re*c[2][0].re + b[0][0].im*c[2][0].im \
                                          +b[0][1].re*c[2][1].re + b[0][1].im*c[2][1].im \
                                          +b[0][2].re*c[2][2].re + b[0][2].im*c[2][2].im;\
                             a[0][2].im = -b[0][0].re*c[2][0].im + b[0][0].im*c[2][0].re \
                                          -b[0][1].re*c[2][1].im + b[0][1].im*c[2][1].re \
                                          -b[0][2].re*c[2][2].im + b[0][2].im*c[2][2].re;\
                                                                                         \
                             a[1][0].re =  b[1][0].re*c[0][0].re + b[1][0].im*c[0][0].im \
                                          +b[1][1].re*c[0][1].re + b[1][1].im*c[0][1].im \
                                          +b[1][2].re*c[0][2].re + b[1][2].im*c[0][2].im;\
                             a[1][0].im = -b[1][0].re*c[0][0].im + b[1][0].im*c[0][0].re \
                                          -b[1][1].re*c[0][1].im + b[1][1].im*c[0][1].re \
                                          -b[1][2].re*c[0][2].im + b[1][2].im*c[0][2].re;\
                             a[1][1].re =  b[1][0].re*c[1][0].re + b[1][0].im*c[1][0].im \
                                          +b[1][1].re*c[1][1].re + b[1][1].im*c[1][1].im \
                                          +b[1][2].re*c[1][2].re + b[1][2].im*c[1][2].im;\
                             a[1][1].im = -b[1][0].re*c[1][0].im + b[1][0].im*c[1][0].re \
                                          -b[1][1].re*c[1][1].im + b[1][1].im*c[1][1].re \
                                          -b[1][2].re*c[1][2].im + b[1][2].im*c[1][2].re;\
                             a[1][2].re =  b[1][0].re*c[2][0].re + b[1][0].im*c[2][0].im \
                                          +b[1][1].re*c[2][1].re + b[1][1].im*c[2][1].im \
                                          +b[1][2].re*c[2][2].re + b[1][2].im*c[2][2].im;\
                             a[1][2].im = -b[1][0].re*c[2][0].im + b[1][0].im*c[2][0].re \
                                          -b[1][1].re*c[2][1].im + b[1][1].im*c[2][1].re \
                                          -b[1][2].re*c[2][2].im + b[1][2].im*c[2][2].re;\
                                                                                         \
                             a[2][0].re =  b[2][0].re*c[0][0].re + b[2][0].im*c[0][0].im \
                                          +b[2][1].re*c[0][1].re + b[2][1].im*c[0][1].im \
                                          +b[2][2].re*c[0][2].re + b[2][2].im*c[0][2].im;\
                             a[2][0].im = -b[2][0].re*c[0][0].im + b[2][0].im*c[0][0].re \
                                          -b[2][1].re*c[0][1].im + b[2][1].im*c[0][1].re \
                                          -b[2][2].re*c[0][2].im + b[2][2].im*c[0][2].re;\
                             a[2][1].re =  b[2][0].re*c[1][0].re + b[2][0].im*c[1][0].im \
                                          +b[2][1].re*c[1][1].re + b[2][1].im*c[1][1].im \
                                          +b[2][2].re*c[1][2].re + b[2][2].im*c[1][2].im;\
                             a[2][1].im = -b[2][0].re*c[1][0].im + b[2][0].im*c[1][0].re \
                                          -b[2][1].re*c[1][1].im + b[2][1].im*c[1][1].re \
                                          -b[2][2].re*c[1][2].im + b[2][2].im*c[1][2].re;\
                             a[2][2].re =  b[2][0].re*c[2][0].re + b[2][0].im*c[2][0].im \
                                          +b[2][1].re*c[2][1].re + b[2][1].im*c[2][1].im \
                                          +b[2][2].re*c[2][2].re + b[2][2].im*c[2][2].im;\
                             a[2][2].im = -b[2][0].re*c[2][0].im + b[2][0].im*c[2][0].re \
                                          -b[2][1].re*c[2][1].im + b[2][1].im*c[2][1].re \
                                          -b[2][2].re*c[2][2].im + b[2][2].im*c[2][2].re;\
})


/* a = adjoint(b)*c, where a,b,c are 3x3 matrices */ 
#define qcd_ADJOINTMUL3x3(a,b,c) ({\
                             a[0][0].re =  b[0][0].re*c[0][0].re + b[0][0].im*c[0][0].im \
                                          +b[1][0].re*c[1][0].re + b[1][0].im*c[1][0].im \
                                          +b[2][0].re*c[2][0].re + b[2][0].im*c[2][0].im;\
                             a[0][0].im =  b[0][0].re*c[0][0].im - b[0][0].im*c[0][0].re \
                                          +b[1][0].re*c[1][0].im - b[1][0].im*c[1][0].re \
                                          +b[2][0].re*c[2][0].im - b[2][0].im*c[2][0].re;\
                             a[0][1].re =  b[0][0].re*c[0][1].re + b[0][0].im*c[0][1].im \
                                          +b[1][0].re*c[1][1].re + b[1][0].im*c[1][1].im \
                                          +b[2][0].re*c[2][1].re + b[2][0].im*c[2][1].im;\
                             a[0][1].im =  b[0][0].re*c[0][1].im - b[0][0].im*c[0][1].re \
                                          +b[1][0].re*c[1][1].im - b[1][0].im*c[1][1].re \
                                          +b[2][0].re*c[2][1].im - b[2][0].im*c[2][1].re;\
                             a[0][2].re =  b[0][0].re*c[0][2].re + b[0][0].im*c[0][2].im \
                                          +b[1][0].re*c[1][2].re + b[1][0].im*c[1][2].im \
                                          +b[2][0].re*c[2][2].re + b[2][0].im*c[2][2].im;\
                             a[0][2].im =  b[0][0].re*c[0][2].im - b[0][0].im*c[0][2].re \
                                          +b[1][0].re*c[1][2].im - b[1][0].im*c[1][2].re \
                                          +b[2][0].re*c[2][2].im - b[2][0].im*c[2][2].re;\
                                                                                         \
                             a[1][0].re =  b[0][1].re*c[0][0].re + b[0][1].im*c[0][0].im \
                                          +b[1][1].re*c[1][0].re + b[1][1].im*c[1][0].im \
                                          +b[2][1].re*c[2][0].re + b[2][1].im*c[2][0].im;\
                             a[1][0].im =  b[0][1].re*c[0][0].im - b[0][1].im*c[0][0].re \
                                          +b[1][1].re*c[1][0].im - b[1][1].im*c[1][0].re \
                                          +b[2][1].re*c[2][0].im - b[2][1].im*c[2][0].re;\
                             a[1][1].re =  b[0][1].re*c[0][1].re + b[0][1].im*c[0][1].im \
                                          +b[1][1].re*c[1][1].re + b[1][1].im*c[1][1].im \
                                          +b[2][1].re*c[2][1].re + b[2][1].im*c[2][1].im;\
                             a[1][1].im =  b[0][1].re*c[0][1].im - b[0][1].im*c[0][1].re \
                                          +b[1][1].re*c[1][1].im - b[1][1].im*c[1][1].re \
                                          +b[2][1].re*c[2][1].im - b[2][1].im*c[2][1].re;\
                             a[1][2].re =  b[0][1].re*c[0][2].re + b[0][1].im*c[0][2].im \
                                          +b[1][1].re*c[1][2].re + b[1][1].im*c[1][2].im \
                                          +b[2][1].re*c[2][2].re + b[2][1].im*c[2][2].im;\
                             a[1][2].im =  b[0][1].re*c[0][2].im - b[0][1].im*c[0][2].re \
                                          +b[1][1].re*c[1][2].im - b[1][1].im*c[1][2].re \
                                          +b[2][1].re*c[2][2].im - b[2][1].im*c[2][2].re;\
                                                                                         \
                             a[2][0].re =  b[0][2].re*c[0][0].re + b[0][2].im*c[0][0].im \
                                          +b[1][2].re*c[1][0].re + b[1][2].im*c[1][0].im \
                                          +b[2][2].re*c[2][0].re + b[2][2].im*c[2][0].im;\
                             a[2][0].im =  b[0][2].re*c[0][0].im - b[0][2].im*c[0][0].re \
                                          +b[1][2].re*c[1][0].im - b[1][2].im*c[1][0].re \
                                          +b[2][2].re*c[2][0].im - b[2][2].im*c[2][0].re;\
                             a[2][1].re =  b[0][2].re*c[0][1].re + b[0][2].im*c[0][1].im \
                                          +b[1][2].re*c[1][1].re + b[1][2].im*c[1][1].im \
                                          +b[2][2].re*c[2][1].re + b[2][2].im*c[2][1].im;\
                             a[2][1].im =  b[0][2].re*c[0][1].im - b[0][2].im*c[0][1].re \
                                          +b[1][2].re*c[1][1].im - b[1][2].im*c[1][1].re \
                                          +b[2][2].re*c[2][1].im - b[2][2].im*c[2][1].re;\
                             a[2][2].re =  b[0][2].re*c[0][2].re + b[0][2].im*c[0][2].im \
                                          +b[1][2].re*c[1][2].re + b[1][2].im*c[1][2].im \
                                          +b[2][2].re*c[2][2].re + b[2][2].im*c[2][2].im;\
                             a[2][2].im =  b[0][2].re*c[0][2].im - b[0][2].im*c[0][2].re \
                                          +b[1][2].re*c[1][2].im - b[1][2].im*c[1][2].re \
                                          +b[2][2].re*c[2][2].im - b[2][2].im*c[2][2].re;\
})


#define qcd_SU3TRACER(u) ( (qcd_real_8) (u[0][0].re + u[1][1].re + u[2][2].re))

 
/* prototypes */ 

void qcd_zero3x3(qcd_complex_16 a[3][3]);
void qcd_unit3x3(qcd_complex_16 a[3][3]);
void qcd_minusUnit3x3(qcd_complex_16 a[3][3]);
void qcd_sub3x3(qcd_complex_16 difference[3][3], qcd_complex_16 minuend[3][3], qcd_complex_16 subtrahend[3][3]);
void qcd_add3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3]);
void qcd_addAdjoint3x3(qcd_complex_16 sum[3][3], qcd_complex_16 summand1[3][3], qcd_complex_16 summand2[3][3]);
void qcd_mul3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3]);
void qcd_mulAdjoint3x3(qcd_complex_16 product[3][3], qcd_complex_16 factor1[3][3], qcd_complex_16 factor2[3][3]);
void qcd_scale3x3(qcd_complex_16 a[3][3], qcd_real_8 r);
void qcd_cScale3x3(qcd_complex_16 a[3][3], qcd_complex_16 c);
qcd_complex_16 qcd_trace3x3(qcd_complex_16 a[3][3]);
void qcd_copy3x3(qcd_complex_16 dest[3][3], qcd_complex_16 src[3][3]);
void qcd_dagger3x3(qcd_complex_16 a[3][3]);
void qcd_zeroVector(qcd_vector *vec);
void qcd_zeroPropagator(qcd_propagator *prop);
void qcd_zeroGaugeField(qcd_gaugeField *u);
void qcd_setVector(qcd_vector *vec, qcd_complex_16 c);
void qcd_copyVector(qcd_vector *dest, qcd_vector *src);
void qcd_copyGaugeField(qcd_gaugeField *dest, qcd_gaugeField *src);
void qcd_copyPropagatorPropagator(qcd_propagator *prop_new, qcd_propagator *prop_old);
void qcd_copyVectorPropagator(qcd_vector *vec, qcd_propagator *prop, qcd_uint_2 nu, qcd_uint_2 c2);
void qcd_copyPropagatorVector(qcd_propagator *prop, qcd_vector *vec, qcd_uint_2 nu, qcd_uint_2 c2);
void qcd_copyVectorPropagator2(qcd_vector *vec, qcd_propagator *prop, qcd_uint_2 mu, qcd_uint_2 c1);
void qcd_copyPropagatorVector2(qcd_propagator *prop, qcd_vector *vec, qcd_uint_2 nu, qcd_uint_2 c2);
void qcd_mulVectorC(qcd_vector *vec, qcd_complex_16 c);
void qcd_mulVectorC3d(qcd_vector *vec, qcd_complex_16 c, qcd_uint_4 t);
void qcd_mulPropagatorC(qcd_propagator *prop, qcd_complex_16 c);
void qcd_mulPropagatorC3d(qcd_propagator *prop, qcd_complex_16 c, qcd_uint_4 t);
void qcd_scaleVector(qcd_vector *vec, qcd_real_8 r);
void qcd_scaleVector3d(qcd_vector *vec, qcd_real_8 r, qcd_uint_4 tt);
void qcd_subVector(qcd_vector *difference, qcd_vector *minuend, qcd_vector *subtrahend);
void qcd_subGaugeField(qcd_gaugeField *difference, qcd_gaugeField *minuend, qcd_gaugeField *subtrahend);
void qcd_addGaugeField(qcd_gaugeField *sum, qcd_gaugeField *summand1, qcd_gaugeField *summand2);
void qcd_scaleGaugeField(qcd_gaugeField *u, qcd_real_8 alpha);
void qcd_addVector(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2);
void qcd_axpyVector(qcd_vector *axpy, qcd_complex_16 a, qcd_vector *x, qcd_vector *y);
void qcd_addVector3d(qcd_vector *sum, qcd_vector *summand1, qcd_vector *summand2, qcd_uint_4 tt);
qcd_complex_16 qcd_mulAdjointVectorVector(qcd_vector *factor1, qcd_vector *factor2);
qcd_real_8 qcd_normVector(qcd_vector *vec);
qcd_real_8 qcd_normGaugeField(qcd_gaugeField *u);
void qcd_conjPropagator(qcd_propagator *prop);
void qcd_projectSU3(qcd_gaugeField *gf);
void qcd_projectSU33d(qcd_gaugeField *gf);
void qcd_projectSU33x3(qcd_complex_16 g[3][3]);
void qcd_randomGaugeField(qcd_gaugeField *u);

#endif
