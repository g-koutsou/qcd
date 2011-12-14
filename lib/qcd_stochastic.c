#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <qcd.h>

void qcd_z2Vector(qcd_vector *vec, int t_slice)
{
  qcd_complex_16 z2[] = {
    (qcd_complex_16){+1., +0.},                                                  
    (qcd_complex_16){+0., +1.},
    (qcd_complex_16){-1., +0.},
    (qcd_complex_16){+0., -1.}};

  qcd_uint_2 *ldims = vec->geo->lL;

  int proc_t = t_slice / ldims[0];
  int my_t_coord = vec->geo->Pos[0];
  int t_slice_loc = t_slice % ldims[0];
  
  memset(&(vec->D[0][0][0].re), 0, 4*3*(vec->geo)->lV*sizeof(qcd_complex_16));
  if(t_slice >= 0) {
    if(my_t_coord == proc_t) { 
      for(int z=0; z<ldims[3]; z++)
	for(int y=0; y<ldims[2]; y++)
	  for(int x=0; x<ldims[1]; x++)
	    {
	      int iv = qcd_LEXIC(t_slice_loc, x, y, z, ldims);
	      for(int sp=0; sp<4; sp++)
		for(int col=0; col<3; col++)
		  {
		    unsigned long int r = qcd_get_rand_4();
		    vec->D[iv][sp][col] = z2[r];
		  }
	    }
    }
  } else if(t_slice < 0) {
    for(int t=0; t<ldims[0]; t++)
      for(int z=0; z<ldims[3]; z++)
	for(int y=0; y<ldims[2]; y++)
	  for(int x=0; x<ldims[1]; x++)
	    {
	      int iv = qcd_LEXIC(t, x, y, z, ldims);
	      for(int sp=0; sp<4; sp++)
		for(int col=0; col<3; col++)
		  {
		    unsigned long int r = qcd_get_rand_4();
		    vec->D[iv][sp][col] = z2[r];
		  }
	    }
  }
  return;
} //end qcd_z2Vector
