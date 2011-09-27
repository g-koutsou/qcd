/* qcd_wilson.c
 * 
 * collection of Wilson Dirac operators
 *
 * Tomasz Korzec 2008
 **********************************************/
 
/* some macros 
   these are valid in the TM chiral representation
   assuming that x and y are 4*3*2 component vectors
   with index i= (is_imaginary)+2*color+2*3*dirac 
*/

#define qcd_ONE_PLUS_GAMMA0(x,y) \
  ({ x[0] = y[0] - y[12]; \
     x[1] = y[1] - y[13]; \
     x[2] = y[2] - y[14]; \
     x[3] = y[3] - y[15]; \
     x[4] = y[4] - y[16]; \
     x[5] = y[5] - y[17]; \
                          \
     x[6] = y[6] - y[18]; \
     x[7] = y[7] - y[19]; \
     x[8] = y[8] - y[20]; \
     x[9] = y[9] - y[21]; \
     x[10]= y[10]- y[22]; \
     x[11]= y[11]- y[23]; \
                          \
     x[12]= -x[0]; \
     x[13]= -x[1]; \
     x[14]= -x[2]; \
     x[15]= -x[3]; \
     x[16]= -x[4]; \
     x[17]= -x[5]; \
                   \
     x[18]= -x[6]; \
     x[19]= -x[7]; \
     x[20]= -x[8]; \
     x[21]= -x[9]; \
     x[22]= -x[10];\
     x[23]= -x[11];\
   })   
   
#define qcd_ONE_PLUS_GAMMA1(x,y) \
  ({ x[0] = y[0] + y[19]; \
     x[1] = y[1] - y[18]; \
     x[2] = y[2] + y[21]; \
     x[3] = y[3] - y[20]; \
     x[4] = y[4] + y[23]; \
     x[5] = y[5] - y[22]; \
                          \
     x[6] = y[6] + y[13]; \
     x[7] = y[7] - y[12]; \
     x[8] = y[8] + y[15]; \
     x[9] = y[9] - y[14]; \
     x[10]= y[10]+ y[17]; \
     x[11]= y[11]- y[16]; \
                          \
     x[12]= -x[7]; \
     x[13]=  x[6]; \
     x[14]= -x[9]; \
     x[15]=  x[8]; \
     x[16]= -x[11];\
     x[17]=  x[10];\
                   \
     x[18]= -x[1]; \
     x[19]=  x[0]; \
     x[20]= -x[3]; \
     x[21]=  x[2]; \
     x[22]= -x[5]; \
     x[23]=  x[4]; \
   })
   
#define qcd_ONE_PLUS_GAMMA2(x,y) \
  ({ x[0] = y[0] - y[18]; \
     x[1] = y[1] - y[19]; \
     x[2] = y[2] - y[20]; \
     x[3] = y[3] - y[21]; \
     x[4] = y[4] - y[22]; \
     x[5] = y[5] - y[23]; \
                          \
     x[6] = y[6] + y[12]; \
     x[7] = y[7] + y[13]; \
     x[8] = y[8] + y[14]; \
     x[9] = y[9] + y[15]; \
     x[10]= y[10]+ y[16]; \
     x[11]= y[11]+ y[17]; \
                          \
     x[12]=   x[6]; \
     x[13]=   x[7]; \
     x[14]=   x[8]; \
     x[15]=   x[9]; \
     x[16]=   x[10];\
     x[17]=   x[11];\
                    \
     x[18]=  -x[0]; \
     x[19]=  -x[1]; \
     x[20]=  -x[2]; \
     x[21]=  -x[3]; \
     x[22]=  -x[4]; \
     x[23]=  -x[5]; \
   }) 
   
#define qcd_ONE_PLUS_GAMMA3(x,y) \
  ({ x[0] = y[0] + y[13]; \
     x[1] = y[1] - y[12]; \
     x[2] = y[2] + y[15]; \
     x[3] = y[3] - y[14]; \
     x[4] = y[4] + y[17]; \
     x[5] = y[5] - y[16]; \
                          \
     x[6] = y[6] - y[19]; \
     x[7] = y[7] + y[18]; \
     x[8] = y[8] - y[21]; \
     x[9] = y[9] + y[20]; \
     x[10]= y[10]- y[23]; \
     x[11]= y[11]+ y[22]; \
                          \
     x[12]=  -x[1]; \
     x[13]=   x[0]; \
     x[14]=  -x[3]; \
     x[15]=   x[2]; \
     x[16]=  -x[5]; \
     x[17]=   x[4]; \
                    \
     x[18]=   x[7]; \
     x[19]=  -x[6]; \
     x[20]=   x[9]; \
     x[21]=  -x[8]; \
     x[22]=   x[11];\
     x[23]=  -x[10];\
   })             
 

#define qcd_ONE_MINUS_GAMMA0(x,y) \
  ({ x[0] = y[0] + y[12]; \
     x[1] = y[1] + y[13]; \
     x[2] = y[2] + y[14]; \
     x[3] = y[3] + y[15]; \
     x[4] = y[4] + y[16]; \
     x[5] = y[5] + y[17]; \
                          \
     x[6] = y[6] + y[18]; \
     x[7] = y[7] + y[19]; \
     x[8] = y[8] + y[20]; \
     x[9] = y[9] + y[21]; \
     x[10]= y[10]+ y[22]; \
     x[11]= y[11]+ y[23]; \
                          \
     x[12]=  x[0]; \
     x[13]=  x[1]; \
     x[14]=  x[2]; \
     x[15]=  x[3]; \
     x[16]=  x[4]; \
     x[17]=  x[5]; \
                   \
     x[18]=  x[6]; \
     x[19]=  x[7]; \
     x[20]=  x[8]; \
     x[21]=  x[9]; \
     x[22]=  x[10];\
     x[23]=  x[11];\
   }) 
 
#define qcd_ONE_MINUS_GAMMA1(x,y) \
  ({ x[0] = y[0] - y[19]; \
     x[1] = y[1] + y[18]; \
     x[2] = y[2] - y[21]; \
     x[3] = y[3] + y[20]; \
     x[4] = y[4] - y[23]; \
     x[5] = y[5] + y[22]; \
                          \
     x[6] = y[6] - y[13]; \
     x[7] = y[7] + y[12]; \
     x[8] = y[8] - y[15]; \
     x[9] = y[9] + y[14]; \
     x[10]= y[10]- y[17]; \
     x[11]= y[11]+ y[16]; \
                          \
     x[12]=  x[7]; \
     x[13]= -x[6]; \
     x[14]=  x[9]; \
     x[15]= -x[8]; \
     x[16]=  x[11];\
     x[17]= -x[10];\
                   \
     x[18]=  x[1]; \
     x[19]= -x[0]; \
     x[20]=  x[3]; \
     x[21]= -x[2]; \
     x[22]=  x[5]; \
     x[23]= -x[4]; \
   }) 
 

#define qcd_ONE_MINUS_GAMMA2(x,y) \
  ({ x[0] = y[0] + y[18]; \
     x[1] = y[1] + y[19]; \
     x[2] = y[2] + y[20]; \
     x[3] = y[3] + y[21]; \
     x[4] = y[4] + y[22]; \
     x[5] = y[5] + y[23]; \
                          \
     x[6] = y[6] - y[12]; \
     x[7] = y[7] - y[13]; \
     x[8] = y[8] - y[14]; \
     x[9] = y[9] - y[15]; \
     x[10]= y[10]- y[16]; \
     x[11]= y[11]- y[17]; \
                          \
     x[12]=  -x[6]; \
     x[13]=  -x[7]; \
     x[14]=  -x[8]; \
     x[15]=  -x[9]; \
     x[16]=  -x[10];\
     x[17]=  -x[11];\
                   \
     x[18]=  x[0]; \
     x[19]=  x[1]; \
     x[20]=  x[2]; \
     x[21]=  x[3]; \
     x[22]=  x[4]; \
     x[23]=  x[5]; \
   })
   
#define qcd_ONE_MINUS_GAMMA3(x,y) \
  ({ x[0] = y[0] - y[13]; \
     x[1] = y[1] + y[12]; \
     x[2] = y[2] - y[15]; \
     x[3] = y[3] + y[14]; \
     x[4] = y[4] - y[17]; \
     x[5] = y[5] + y[16]; \
                          \
     x[6] = y[6] + y[19]; \
     x[7] = y[7] - y[18]; \
     x[8] = y[8] + y[21]; \
     x[9] = y[9] - y[20]; \
     x[10]= y[10]+ y[23]; \
     x[11]= y[11]- y[22]; \
                          \
     x[12]=   x[1]; \
     x[13]=  -x[0]; \
     x[14]=   x[3]; \
     x[15]=  -x[2]; \
     x[16]=   x[5]; \
     x[17]=  -x[4]; \
                    \
     x[18]=  -x[7]; \
     x[19]=   x[6]; \
     x[20]=  -x[9]; \
     x[21]=   x[8]; \
     x[22]=  -x[11];\
     x[23]=   x[10];\
   })
   
#define qcd_MASSTERMS(x,y,mp4,mu) \
  ({ x[0] = mp4*y[0]-mu*y[1];  \
     x[1] = mp4*y[1]+mu*y[0];  \
     x[2] = mp4*y[2]-mu*y[3];  \
     x[3] = mp4*y[3]+mu*y[2];  \
     x[4] = mp4*y[4]-mu*y[5];  \
     x[5] = mp4*y[5]+mu*y[4];  \
                               \
     x[6] = mp4*y[6] -mu*y[7]; \
     x[7] = mp4*y[7] +mu*y[6]; \
     x[8] = mp4*y[8] -mu*y[9]; \
     x[9] = mp4*y[9] +mu*y[8]; \
     x[10]= mp4*y[10]-mu*y[11];\
     x[11]= mp4*y[11]+mu*y[10];\
                               \
     x[12]= mp4*y[12]+mu*y[13];\
     x[13]= mp4*y[13]-mu*y[12];\
     x[14]= mp4*y[14]+mu*y[15];\
     x[15]= mp4*y[15]-mu*y[14];\
     x[16]= mp4*y[16]+mu*y[17];\
     x[17]= mp4*y[17]-mu*y[16];\
                               \
     x[18]= mp4*y[18]+mu*y[19];\
     x[19]= mp4*y[19]-mu*y[18];\
     x[20]= mp4*y[20]+mu*y[21];\
     x[21]= mp4*y[21]-mu*y[20];\
     x[22]= mp4*y[22]+mu*y[23];\
     x[23]= mp4*y[23]-mu*y[22];\
   }) 
   
#define qcd_APPLY_U(u,x,y)  \
   ({ x[0] = u[0] *y[0]  - u[1] *y[1]  + u[2] *y[2]  - u[3] *y[3]  + u[4] *y[4]  - u[5] *y[5];     \
      x[1] = u[0] *y[1]  + u[1] *y[0]  + u[2] *y[3]  + u[3] *y[2]  + u[4] *y[5]  + u[5] *y[4];     \
      x[2] = u[6] *y[0]  - u[7] *y[1]  + u[8] *y[2]  - u[9] *y[3]  + u[10]*y[4]  - u[11]*y[5];     \
      x[3] = u[6] *y[1]  + u[7] *y[0]  + u[8] *y[3]  + u[9] *y[2]  + u[10]*y[5]  + u[11]*y[4];     \
      x[4] = u[12]*y[0]  - u[13]*y[1]  + u[14]*y[2]  - u[15]*y[3]  + u[16]*y[4]  - u[17]*y[5];     \
      x[5] = u[12]*y[1]  + u[13]*y[0]  + u[14]*y[3]  + u[15]*y[2]  + u[16]*y[5]  + u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  - u[1] *y[7]  + u[2] *y[8]  - u[3] *y[9]  + u[4] *y[10] - u[5] *y[11];    \
      x[7] = u[0] *y[7]  + u[1] *y[6]  + u[2] *y[9]  + u[3] *y[8]  + u[4] *y[11] + u[5] *y[10];    \
      x[8] = u[6] *y[6]  - u[7] *y[7]  + u[8] *y[8]  - u[9] *y[9]  + u[10]*y[10] - u[11]*y[11];    \
      x[9] = u[6] *y[7]  + u[7] *y[6]  + u[8] *y[9]  + u[9] *y[8]  + u[10]*y[11] + u[11]*y[10];    \
      x[10]= u[12]*y[6]  - u[13]*y[7]  + u[14]*y[8]  - u[15]*y[9]  + u[16]*y[10] - u[17]*y[11];    \
      x[11]= u[12]*y[7]  + u[13]*y[6]  + u[14]*y[9]  + u[15]*y[8]  + u[16]*y[11] + u[17]*y[10];    \
                                                                                                   \
      x[12]= u[0] *y[12] - u[1] *y[13] + u[2] *y[14] - u[3] *y[15] + u[4] *y[16] - u[5] *y[17];    \
      x[13]= u[0] *y[13] + u[1] *y[12] + u[2] *y[15] + u[3] *y[14] + u[4] *y[17] + u[5] *y[16];    \
      x[14]= u[6] *y[12] - u[7] *y[13] + u[8] *y[14] - u[9] *y[15] + u[10]*y[16] - u[11]*y[17];    \
      x[15]= u[6] *y[13] + u[7] *y[12] + u[8] *y[15] + u[9] *y[14] + u[10]*y[17] + u[11]*y[16];    \
      x[16]= u[12]*y[12] - u[13]*y[13] + u[14]*y[14] - u[15]*y[15] + u[16]*y[16] - u[17]*y[17];    \
      x[17]= u[12]*y[13] + u[13]*y[12] + u[14]*y[15] + u[15]*y[14] + u[16]*y[17] + u[17]*y[16];    \
                                                                                                   \
      x[18]= u[0] *y[18] - u[1] *y[19] + u[2] *y[20] - u[3] *y[21] + u[4] *y[22] - u[5] *y[23];    \
      x[19]= u[0] *y[19] + u[1] *y[18] + u[2] *y[21] + u[3] *y[20] + u[4] *y[23] + u[5] *y[22];    \
      x[20]= u[6] *y[18] - u[7] *y[19] + u[8] *y[20] - u[9] *y[21] + u[10]*y[22] - u[11]*y[23];    \
      x[21]= u[6] *y[19] + u[7] *y[18] + u[8] *y[21] + u[9] *y[20] + u[10]*y[23] + u[11]*y[22];    \
      x[22]= u[12]*y[18] - u[13]*y[19] + u[14]*y[20] - u[15]*y[21] + u[16]*y[22] - u[17]*y[23];    \
      x[23]= u[12]*y[19] + u[13]*y[18] + u[14]*y[21] + u[15]*y[20] + u[16]*y[23] + u[17]*y[22];    \
   })
   
#define qcd_APPLY_U_DAGGER(u,x,y)  \
   ({ x[0] = u[0] *y[0]  + u[1] *y[1]  + u[6] *y[2]  + u[7] *y[3]  + u[12]*y[4]  + u[13]*y[5];     \
      x[1] = u[0] *y[1]  - u[1] *y[0]  + u[6] *y[3]  - u[7] *y[2]  + u[12]*y[5]  - u[13]*y[4];     \
      x[2] = u[2] *y[0]  + u[3] *y[1]  + u[8] *y[2]  + u[9] *y[3]  + u[14]*y[4]  + u[15]*y[5];     \
      x[3] = u[2] *y[1]  - u[3] *y[0]  + u[8] *y[3]  - u[9] *y[2]  + u[14]*y[5]  - u[15]*y[4];     \
      x[4] = u[4] *y[0]  + u[5] *y[1]  + u[10]*y[2]  + u[11]*y[3]  + u[16]*y[4]  + u[17]*y[5];     \
      x[5] = u[4] *y[1]  - u[5] *y[0]  + u[10]*y[3]  - u[11]*y[2]  + u[16]*y[5]  - u[17]*y[4];     \
                                                                                                   \
      x[6] = u[0] *y[6]  + u[1] *y[7]  + u[6] *y[8]  + u[7] *y[9]  + u[12]*y[10] + u[13]*y[11];    \
      x[7] = u[0] *y[7]  - u[1] *y[6]  + u[6] *y[9]  - u[7] *y[8]  + u[12]*y[11] - u[13]*y[10];    \
      x[8] = u[2] *y[6]  + u[3] *y[7]  + u[8] *y[8]  + u[9] *y[9]  + u[14]*y[10] + u[15]*y[11];    \
      x[9] = u[2] *y[7]  - u[3] *y[6]  + u[8] *y[9]  - u[9] *y[8]  + u[14]*y[11] - u[15]*y[10];    \
      x[10]= u[4] *y[6]  + u[5] *y[7]  + u[10]*y[8]  + u[11]*y[9]  + u[16]*y[10] + u[17]*y[11];    \
      x[11]= u[4] *y[7]  - u[5] *y[6]  + u[10]*y[9]  - u[11]*y[8]  + u[16]*y[11] - u[17]*y[10];    \
                                                                                                   \
      x[12]= u[0] *y[12] + u[1] *y[13] + u[6] *y[14] + u[7] *y[15] + u[12]*y[16] + u[13]*y[17];    \
      x[13]= u[0] *y[13] - u[1] *y[12] + u[6] *y[15] - u[7] *y[14] + u[12]*y[17] - u[13]*y[16];    \
      x[14]= u[2] *y[12] + u[3] *y[13] + u[8] *y[14] + u[9] *y[15] + u[14]*y[16] + u[15]*y[17];    \
      x[15]= u[2] *y[13] - u[3] *y[12] + u[8] *y[15] - u[9] *y[14] + u[14]*y[17] - u[15]*y[16];    \
      x[16]= u[4] *y[12] + u[5] *y[13] + u[10]*y[14] + u[11]*y[15] + u[16]*y[16] + u[17]*y[17];    \
      x[17]= u[4] *y[13] - u[5] *y[12] + u[10]*y[15] - u[11]*y[14] + u[16]*y[17] - u[17]*y[16];    \
                                                                                                   \
      x[18]= u[0] *y[18] + u[1] *y[19] + u[6] *y[20] + u[7] *y[21] + u[12]*y[22] + u[13]*y[23];    \
      x[19]= u[0] *y[19] - u[1] *y[18] + u[6] *y[21] - u[7] *y[20] + u[12]*y[23] - u[13]*y[22];    \
      x[20]= u[2] *y[18] + u[3] *y[19] + u[8] *y[20] + u[9] *y[21] + u[14]*y[22] + u[15]*y[23];    \
      x[21]= u[2] *y[19] - u[3] *y[18] + u[8] *y[21] - u[9] *y[20] + u[14]*y[23] - u[15]*y[22];    \
      x[22]= u[4] *y[18] + u[5] *y[19] + u[10]*y[20] + u[11]*y[21] + u[16]*y[22] + u[17]*y[23];    \
      x[23]= u[4] *y[19] - u[5] *y[18] + u[10]*y[21] - u[11]*y[20] + u[16]*y[23] - u[17]*y[22];    \
   }) 
   
   
#define qcd_SUM_OP(x,a,b,c,d,e,f,g,h,i)  \
   ({ x[0] = a[0] -0.5*(b[0] +c[0] +d[0] +e[0] +f[0] +g[0] +h[0] +i[0]);        \
      x[1] = a[1] -0.5*(b[1] +c[1] +d[1] +e[1] +f[1] +g[1] +h[1] +i[1]);        \
      x[2] = a[2] -0.5*(b[2] +c[2] +d[2] +e[2] +f[2] +g[2] +h[2] +i[2]);        \
      x[3] = a[3] -0.5*(b[3] +c[3] +d[3] +e[3] +f[3] +g[3] +h[3] +i[3]);        \
      x[4] = a[4] -0.5*(b[4] +c[4] +d[4] +e[4] +f[4] +g[4] +h[4] +i[4]);        \
      x[5] = a[5] -0.5*(b[5] +c[5] +d[5] +e[5] +f[5] +g[5] +h[5] +i[5]);        \
      x[6] = a[6] -0.5*(b[6] +c[6] +d[6] +e[6] +f[6] +g[6] +h[6] +i[6]);        \
      x[7] = a[7] -0.5*(b[7] +c[7] +d[7] +e[7] +f[7] +g[7] +h[7] +i[7]);        \
      x[8] = a[8] -0.5*(b[8] +c[8] +d[8] +e[8] +f[8] +g[8] +h[8] +i[8]);        \
      x[9] = a[9] -0.5*(b[9] +c[9] +d[9] +e[9] +f[9] +g[9] +h[9] +i[9]);        \
      x[10]= a[10]-0.5*(b[10]+c[10]+d[10]+e[10]+f[10]+g[10]+h[10]+i[10]);       \
      x[11]= a[11]-0.5*(b[11]+c[11]+d[11]+e[11]+f[11]+g[11]+h[11]+i[11]);       \
      x[12]= a[12]-0.5*(b[12]+c[12]+d[12]+e[12]+f[12]+g[12]+h[12]+i[12]);       \
      x[13]= a[13]-0.5*(b[13]+c[13]+d[13]+e[13]+f[13]+g[13]+h[13]+i[13]);       \
      x[14]= a[14]-0.5*(b[14]+c[14]+d[14]+e[14]+f[14]+g[14]+h[14]+i[14]);       \
      x[15]= a[15]-0.5*(b[15]+c[15]+d[15]+e[15]+f[15]+g[15]+h[15]+i[15]);       \
      x[16]= a[16]-0.5*(b[16]+c[16]+d[16]+e[16]+f[16]+g[16]+h[16]+i[16]);       \
      x[17]= a[17]-0.5*(b[17]+c[17]+d[17]+e[17]+f[17]+g[17]+h[17]+i[17]);       \
      x[18]= a[18]-0.5*(b[18]+c[18]+d[18]+e[18]+f[18]+g[18]+h[18]+i[18]);       \
      x[19]= a[19]-0.5*(b[19]+c[19]+d[19]+e[19]+f[19]+g[19]+h[19]+i[19]);       \
      x[20]= a[20]-0.5*(b[20]+c[20]+d[20]+e[20]+f[20]+g[20]+h[20]+i[20]);       \
      x[21]= a[21]-0.5*(b[21]+c[21]+d[21]+e[21]+f[21]+g[21]+h[21]+i[21]);       \
      x[22]= a[22]-0.5*(b[22]+c[22]+d[22]+e[22]+f[22]+g[22]+h[22]+i[22]);       \
      x[23]= a[23]-0.5*(b[23]+c[23]+d[23]+e[23]+f[23]+g[23]+h[23]+i[23]);       \
   })
   
#define qcd_PHASE_MUL(x,t,pr,pi)   \
   ({ t    =x[0] *pr -x[1] *pi;    \
      x[1] =x[0] *pi +x[1] *pr;    \
      x[0] =t;                     \
      t    =x[2] *pr -x[3] *pi;    \
      x[3] =x[2] *pi +x[3] *pr;    \
      x[2] =t;                     \
      t    =x[4] *pr -x[5] *pi;    \
      x[5] =x[4] *pi +x[5] *pr;    \
      x[4] =t;                     \
      t    =x[6] *pr -x[7] *pi;    \
      x[7] =x[6] *pi +x[7] *pr;    \
      x[6] =t;                     \
      t    =x[8] *pr -x[9] *pi;    \
      x[9] =x[8] *pi +x[9] *pr;    \
      x[8] =t;                     \
      t    =x[10]*pr -x[11]*pi;    \
      x[11]=x[10]*pi +x[11]*pr;    \
      x[10]=t;                     \
      t    =x[12]*pr -x[13]*pi;    \
      x[13]=x[12]*pi +x[13]*pr;    \
      x[12]=t;                     \
      t    =x[14]*pr -x[15]*pi;    \
      x[15]=x[14]*pi +x[15]*pr;    \
      x[14]=t;                     \
      t    =x[16]*pr -x[17]*pi;    \
      x[17]=x[16]*pi +x[17]*pr;    \
      x[16]=t;                     \
      t    =x[18]*pr -x[19]*pi;    \
      x[19]=x[18]*pi +x[19]*pr;    \
      x[18]=t;                     \
      t    =x[20]*pr -x[21]*pi;    \
      x[21]=x[20]*pi +x[21]*pr;    \
      x[20]=t;                     \
      t    =x[22]*pr -x[23]*pi;    \
      x[23]=x[22]*pi +x[23]*pr;    \
      x[22]=t;                     \
   })
   
   
/* apply the Wilson TM dirac operator on vector v 
 * u  gauge field
 * m  mass parameter
 * tm twisted mass parameter 
 */
int qcd_applyWilsonTMOp(qcd_vector *Dv, qcd_vector *v, qcd_gaugeField *u, qcd_real_8 m, qcd_real_8 tm)
{
   qcd_uint_2 b;

   qcd_uint_4 t,x,y,z;
   qcd_uint_8 l,j;
   qcd_uint_4 b0,b1,b2,b3;
   qcd_uint_2 bl0[4] = {1,0,0,0}; 
   qcd_uint_2 bl1[4] = {2,2,1,1};
   qcd_uint_2 bl2[4] = {3,3,3,2};
   qcd_uint_2 isPeriodic[4];
   
   qcd_real_8 epg0[24], emg0[24], epg1[24], emg1[24]; // vectors to store (1+-gamma_mu)psi
   qcd_real_8 epg2[24], emg2[24], epg3[24], emg3[24];
   qcd_real_8 uepg0[24], uemg0[24], uepg1[24], uemg1[24];
   qcd_real_8 uepg2[24], uemg2[24], uepg3[24], uemg3[24];
   qcd_real_8 diag[24];
   qcd_real_8 *psi;
   qcd_real_8 mplus4=4.0+m;
   qcd_real_8 *uu;
   qcd_real_8 tmp;
   
   qcd_complex_16 pphases[4];
   qcd_complex_16 mphases[4];
   
   if(!v->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Vector must be properly initialized\n");
      return(1);
   }
   if(!u->initialized)
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Gauge field must be properly initialized\n");
      return(1);
   }
   if(v->geo->initialized!=1 ||
      u->geo->initialized!=1 ||
      v->geo->L[0] != u->geo->L[0] ||
      v->geo->L[1] != u->geo->L[1] ||
      v->geo->L[2] != u->geo->L[2] ||
      v->geo->L[3] != u->geo->L[3] ||
      v->geo->lL[0] != u->geo->lL[0] ||
      v->geo->lL[1] != u->geo->lL[1] ||
      v->geo->lL[2] != u->geo->lL[2] ||
      v->geo->lL[3] != u->geo->lL[3] ||
      v->geo->Pos[0] != u->geo->Pos[0] ||
      v->geo->Pos[1] != u->geo->Pos[1] ||
      v->geo->Pos[2] != u->geo->Pos[2] ||
      v->geo->Pos[3] != u->geo->Pos[3])
   {
      fprintf(stderr,"Error in qcd_ApplyWilsonTMOp! Geometry not initialized or missmatched\n");
      return(1);   
   }       
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication\n");
   qcd_communicateVectorPMGaugeP(v,u);   
      
      
   for(b=0; b<4; b++)
   {
      if(u->geo->theta[b] != 0)
      {
         pphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/u->geo->L[b]),  sin(u->geo->theta[b]/u->geo->L[b])};
         mphases[b] = (qcd_complex_16) {cos(u->geo->theta[b]/u->geo->L[b]), -sin(u->geo->theta[b]/u->geo->L[b])};
         isPeriodic[b] = 0;
      }
      else
      {
         pphases[b] = (qcd_complex_16) {1,0};
         mphases[b] = (qcd_complex_16) {1,0};
         isPeriodic[b] = 1;
      }   
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting inner-points loop\n");
   //calculate interior points
   if((u->geo->lL[0] > 2) && (u->geo->lL[1] > 2) && (u->geo->lL[2] > 2) && (u->geo->lL[3] > 2))
   {
      for(z=1; z<u->geo->lL[3]-1; z++)
      for(y=1; y<u->geo->lL[2]-1; y++)
      for(x=1; x<u->geo->lL[1]-1; x++)
      for(t=1; t<u->geo->lL[0]-1; t++)
      {
         l = qcd_LEXIC(t,x,y,z,u->geo->lL);         
         
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][0]][0][0].re);         
         qcd_ONE_PLUS_GAMMA0(epg0,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][1]][0][0].re);         
         qcd_ONE_PLUS_GAMMA1(epg1,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][2]][0][0].re);         
         qcd_ONE_PLUS_GAMMA2(epg2,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->minus[l][3]][0][0].re);         
         qcd_ONE_PLUS_GAMMA3(epg3,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][0]][0][0].re);         
         qcd_ONE_MINUS_GAMMA0(emg0,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][1]][0][0].re);         
         qcd_ONE_MINUS_GAMMA1(emg1,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][2]][0][0].re);         
         qcd_ONE_MINUS_GAMMA2(emg2,psi);
         psi = (qcd_real_8*) &(v->D[v->geo->plus[l][3]][0][0].re);         
         qcd_ONE_MINUS_GAMMA3(emg3,psi); 

         uu = (qcd_real_8*) &(u->D[l][0][0][0].re);
         qcd_APPLY_U(uu,uemg0,emg0);
         uu = (qcd_real_8*) &(u->D[l][1][0][0].re);
         qcd_APPLY_U(uu,uemg1,emg1);
         uu = (qcd_real_8*) &(u->D[l][2][0][0].re);
         qcd_APPLY_U(uu,uemg2,emg2);
         uu = (qcd_real_8*) &(u->D[l][3][0][0].re);
         qcd_APPLY_U(uu,uemg3,emg3);

         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][0]][0][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg0,epg0);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][1]][1][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg1,epg1);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][2]][2][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg2,epg2);
         uu = (qcd_real_8*) &(u->D[u->geo->minus[l][3]][3][0][0].re);
         qcd_APPLY_U_DAGGER(uu,uepg3,epg3);                           
         
         psi = (qcd_real_8*) &(v->D[l][0][0].re);
         qcd_MASSTERMS(diag,psi,mplus4,tm);
         
         psi = (qcd_real_8*) &(Dv->D[l][0][0].re);
         if(!isPeriodic[0])
         {         
            qcd_PHASE_MUL(uepg0,tmp,mphases[0].re,mphases[0].im);
            qcd_PHASE_MUL(uemg0,tmp,pphases[0].re,pphases[0].im); 
         }
         if(!isPeriodic[1])
         {         
            qcd_PHASE_MUL(uepg1,tmp,mphases[1].re,mphases[1].im);
            qcd_PHASE_MUL(uemg1,tmp,pphases[1].re,pphases[1].im); 
         }
         if(!isPeriodic[2])
         {         
            qcd_PHASE_MUL(uepg2,tmp,mphases[2].re,mphases[2].im);
            qcd_PHASE_MUL(uemg2,tmp,pphases[2].re,pphases[2].im); 
         }
         if(!isPeriodic[3])
         {         
            qcd_PHASE_MUL(uepg3,tmp,mphases[3].re,mphases[3].im);
            qcd_PHASE_MUL(uemg3,tmp,pphases[3].re,pphases[3].im); 
         }
         qcd_SUM_OP(psi,diag,uepg0,uemg0,uepg1,uemg1,uepg2,uemg2,uepg3,uemg3);
      }
   }
   
   //if(u->geo->myid==0) printf("Wilson operator: starting communication finalization\n");
   qcd_waitall(v->geo);
   
   
   //if(u->geo->myid==0) printf("Wilson operator: starting boundaries loop\n");
   //calculate boundary points
   for(j=0; j<v->geo->edgePoints; j++)
   {
      l=v->geo->edge[j];  
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][0]][0][0].re);         
      qcd_ONE_PLUS_GAMMA0(epg0,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][1]][0][0].re);         
      qcd_ONE_PLUS_GAMMA1(epg1,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][2]][0][0].re);         
      qcd_ONE_PLUS_GAMMA2(epg2,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->minus[l][3]][0][0].re);         
      qcd_ONE_PLUS_GAMMA3(epg3,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][0]][0][0].re);         
      qcd_ONE_MINUS_GAMMA0(emg0,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][1]][0][0].re);         
      qcd_ONE_MINUS_GAMMA1(emg1,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][2]][0][0].re);         
      qcd_ONE_MINUS_GAMMA2(emg2,psi);
      psi = (qcd_real_8*) &(v->D[v->geo->plus[l][3]][0][0].re);         
      qcd_ONE_MINUS_GAMMA3(emg3,psi); 

      uu = (qcd_real_8*) &(u->D[l][0][0][0].re);
      qcd_APPLY_U(uu,uemg0,emg0);
      uu = (qcd_real_8*) &(u->D[l][1][0][0].re);
      qcd_APPLY_U(uu,uemg1,emg1);
      uu = (qcd_real_8*) &(u->D[l][2][0][0].re);
      qcd_APPLY_U(uu,uemg2,emg2);
      uu = (qcd_real_8*) &(u->D[l][3][0][0].re);
      qcd_APPLY_U(uu,uemg3,emg3);

      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][0]][0][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg0,epg0);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][1]][1][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg1,epg1);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][2]][2][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg2,epg2);
      uu = (qcd_real_8*) &(u->D[u->geo->minus[l][3]][3][0][0].re);
      qcd_APPLY_U_DAGGER(uu,uepg3,epg3);                           

      psi = (qcd_real_8*) &(v->D[l][0][0].re);
      qcd_MASSTERMS(diag,psi,mplus4,tm);

      psi = (qcd_real_8*) &(Dv->D[l][0][0].re);
      if(!isPeriodic[0])
      {         
         qcd_PHASE_MUL(uepg0,tmp,mphases[0].re,mphases[0].im);
         qcd_PHASE_MUL(uemg0,tmp,pphases[0].re,pphases[0].im); 
      }
      if(!isPeriodic[1])
      {         
         qcd_PHASE_MUL(uepg1,tmp,mphases[1].re,mphases[1].im);
         qcd_PHASE_MUL(uemg1,tmp,pphases[1].re,pphases[1].im); 
      }
      if(!isPeriodic[2])
      {         
         qcd_PHASE_MUL(uepg2,tmp,mphases[2].re,mphases[2].im);
         qcd_PHASE_MUL(uemg2,tmp,pphases[2].re,pphases[2].im); 
      }
      if(!isPeriodic[3])
      {         
         qcd_PHASE_MUL(uepg3,tmp,mphases[3].re,mphases[3].im);
         qcd_PHASE_MUL(uemg3,tmp,pphases[3].re,pphases[3].im); 
      }
      qcd_SUM_OP(psi,diag,uepg0,uemg0,uepg1,uemg1,uepg2,uemg2,uepg3,uemg3);                 
   }//end loop over boundary points
   
   
   
   return(0);
}//end applyWilsonTMOp 
