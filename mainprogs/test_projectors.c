#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>
#include <omp.h>
#include <projectors_pb.h>


int main(int argc, char *argv[]){
  
  int i,j,k,i1,i2,i3,al,be,ga,de,prid,projlist[3] = {16,15,4},proj;

  qcd_complex_16 Pr32[9][4][4],delta,Fin_proj[3][4][4][3][3],Wrg_proj[3][4][4][3][3];
  qcd_real_8 fac = 1.0/3.0;
  qcd_complex_16 cfac;
  qcd_complex_16 gamma12[4][4],gamma13[4][4],gamma23[4][4];

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      gamma12[i][j] = (qcd_complex_16) {0,0};
      gamma13[i][j] = (qcd_complex_16) {0,0};
      gamma23[i][j] = (qcd_complex_16) {0,0};
      
      for(k=0;k<4;k++){
	gamma12[i][j] = qcd_CADD(gamma12[i][j],qcd_CMUL(qcd_GAMMA[1][i][k],qcd_GAMMA[2][k][j]));
	gamma13[i][j] = qcd_CADD(gamma13[i][j],qcd_CMUL(qcd_GAMMA[1][i][k],qcd_GAMMA[3][k][j]));
	gamma23[i][j] = qcd_CADD(gamma23[i][j],qcd_CMUL(qcd_GAMMA[2][i][k],qcd_GAMMA[3][k][j]));
      }
    }
  }

  for(i=0;i<9;i++){
    for(al=0;al<4;al++){
      for(be=0;be<4;be++){
        Pr32[i][al][be] = (qcd_complex_16) {0.0,0.0};
      }
    }
  }

  // Define the digonal elements of the projector to 3/2
  // elements of [0][al][be], [4][al][be] [8][al][be] for al!=be are zero by definition 
  delta = (qcd_complex_16) {1.0,0.0};
  cfac  = (qcd_complex_16) {1.0/3.0,0.0};
  for(al=0;al<4;al++){
    Pr32[0][al][al] = qcd_CSUB(delta,cfac); // Pr_11
    Pr32[4][al][al] = qcd_CSUB(delta,cfac); // Pr_22
    Pr32[8][al][al] = qcd_CSUB(delta,cfac); // Pr_33
  }

  // Define the rest elements
  delta = (qcd_complex_16) {0.0,0.0};
  for(al=0;al<4;al++){
    for(be=0;be<4;be++){
      Pr32[1][al][be] = qcd_CSUB(delta,qcd_CSCALE(gamma12[al][be],fac)); // Pr_12
      Pr32[2][al][be] = qcd_CSUB(delta,qcd_CSCALE(gamma13[al][be],fac)); // Pr_13
      Pr32[5][al][be] = qcd_CSUB(delta,qcd_CSCALE(gamma23[al][be],fac)); // Pr_23

      Pr32[3][al][be] = qcd_CSUB(delta,qcd_CSCALE(qcd_CSCALE(gamma12[al][be],-1.0),fac)); // Pr_21
      Pr32[6][al][be] = qcd_CSUB(delta,qcd_CSCALE(qcd_CSCALE(gamma13[al][be],-1.0),fac)); // Pr_31
      Pr32[7][al][be] = qcd_CSUB(delta,qcd_CSCALE(qcd_CSCALE(gamma23[al][be],-1.0),fac)); // Pr_32
    }
  }

  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      for(al=0;al<4;al++){
	for(be=0;be<4;be++){
	  i1 = j+3*i;
	  printf("Pr %d %d %d %d = %lf %lf\n",al+1,be+1,i+1,j+1,Pr32[i1][al][be].re,Pr32[i1][al][be].im);
	}
      }
    }
  }

  printf("\n\n\n");


  //-Define the projector 
  for(prid=0;prid<3;prid++){
    proj = projlist[prid];

    for(al=0;al<4;al++){
      for(de=0;de<4;de++){
	for(k=0;k<3;k++){
	  for(j=0;j<3;j++){
	    Fin_proj[prid][al][de][k][j] = (qcd_complex_16) {0.0,0.0};

	    for(i=0;i<3;i++){
	      i1 = i+3*k;
	      i2 = j+3*i;
	      for(be=0;be<4;be++){
		for(ga=0;ga<4;ga++){

		  Fin_proj[prid][al][de][k][j] = qcd_CADD(Fin_proj[prid][al][de][k][j],qcd_CMUL(Pr32[i1][al][be],qcd_CMUL(PROJECTOR[proj][be][ga],Pr32[i2][ga][de])));

		}//-ga
	      }//-be
	    }//-i
	  }//-j
	}//-k
      }//-de
    }//-al
  }//-prid


  /* for(al=0;al<4;al++){ */
  /*   for(de=0;de<4;de++){ */
  /*     for(k=0;k<3;k++){ */
  /* 	for(j=0;j<3;j++){ */
  /* 	  Fin_proj[0][al][de][k][j] = (qcd_complex_16) {0.0,0.0}; */
  /* 	  Fin_proj[1][al][de][k][j] = (qcd_complex_16) {0.0,0.0}; */
  /* 	  Fin_proj[2][al][de][k][j] = (qcd_complex_16) {0.0,0.0}; */
  /* 	  Wrg_proj[0][al][de][k][j] = (qcd_complex_16) {0.0,0.0}; */
  /* 	  Wrg_proj[1][al][de][k][j] = (qcd_complex_16) {0.0,0.0}; */
  /* 	  Wrg_proj[2][al][de][k][j] = (qcd_complex_16) {0.0,0.0}; */

  /* 	  i3 = j+3*k; */
  /* 	  for(be=0;be<4;be++){ */
  /* 	    Wrg_proj[0][al][de][k][j] = qcd_CADD(Wrg_proj[0][al][de][k][j],qcd_CMUL(Pr32[i3][al][be],PROJECTOR[16][be][de])); */
  /* 	    Wrg_proj[1][al][de][k][j] = qcd_CADD(Wrg_proj[1][al][de][k][j],qcd_CMUL(Pr32[i3][al][be],PROJECTOR[15][be][de])); */
  /* 	    Wrg_proj[2][al][de][k][j] = qcd_CADD(Wrg_proj[2][al][de][k][j],qcd_CMUL(Pr32[i3][al][be],PROJECTOR[4][be][de])); */
  /* 	  } */

  /* 	  for(i=0;i<3;i++){ */
  /* 	    i1 = i+3*k; */
  /* 	    i2 = j+3*i; */
  /* 	    for(be=0;be<4;be++){ */
  /* 	      for(ga=0;ga<4;ga++){ */
  /* 		Fin_proj[0][al][de][k][j] = qcd_CADD(Fin_proj[0][al][de][k][j],qcd_CMUL(Pr32[i1][al][be],qcd_CMUL(PROJECTOR[16][be][ga],Pr32[i2][ga][de]))); */
  /* 		Fin_proj[1][al][de][k][j] = qcd_CADD(Fin_proj[1][al][de][k][j],qcd_CMUL(Pr32[i1][al][be],qcd_CMUL(PROJECTOR[15][be][ga],Pr32[i2][ga][de]))); */
  /* 		Fin_proj[2][al][de][k][j] = qcd_CADD(Fin_proj[2][al][de][k][j],qcd_CMUL(Pr32[i1][al][be],qcd_CMUL(PROJECTOR[4][be][ga],Pr32[i2][ga][de]))); */
  /* 	      } */
  /* 	    } */
  /* 	  } */
  /* 	} */
  /*     } */
  /*   } */
  /* } */


  for(k=0;k<3;k++){
    for(j=0;j<3;j++){
      for(al=0;al<4;al++){
	for(de=0;de<4;de++){
	  printf("(%d,%d,%d,%d) (%lf,%lf) (%lf,%lf) (%lf,%lf)\n",al+1,de+1,k+1,j+1,
		 Fin_proj[0][al][de][k][j].re,Fin_proj[0][al][de][k][j].im,
		 Fin_proj[1][al][de][k][j].re,Fin_proj[1][al][de][k][j].im,
		 Fin_proj[2][al][de][k][j].re,Fin_proj[2][al][de][k][j].im); 
	  /* Wrg_proj[0][al][de][k][j].re,Wrg_proj[0][al][de][k][j].im, */
	  /* Wrg_proj[1][al][de][k][j].re,Wrg_proj[1][al][de][k][j].im, */
	  /* Wrg_proj[2][al][de][k][j].re,Wrg_proj[2][al][de][k][j].im);  */
	}
      }
    }
  }
  

  return 0;
}
