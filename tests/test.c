#include <time.h>
#include <stdio.h>
#include <stdlib.h>

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

void dbyte_swap(R,N)
     int   N;
     char *R;
{
  int  i;
  char swap;

  for(i=0;i<N;i++){
    swap = *(R + i*8    ); *(R + i*8    ) = *(R + i*8 + 7);  *(R + i*8 + 7) = swap;
    swap = *(R + i*8 + 1); *(R + i*8 + 1) = *(R + i*8 + 6);  *(R + i*8 + 6) = swap;
    swap = *(R + i*8 + 2); *(R + i*8 + 2) = *(R + i*8 + 5);  *(R + i*8 + 5) = swap;
    swap = *(R + i*8 + 3); *(R + i*8 + 3) = *(R + i*8 + 4);  *(R + i*8 + 4) = swap;
  };
};

int main()
{
   double *tmp;
   int i;
   clock_t t1,t2;
   double dif;
   tmp = malloc(2500000*sizeof(double));
   for (i=0; i<2500000; i++)
   {
      tmp[i] = i+1/(i+1);
   }
   t1=clock();   
   dbyte_swap(tmp,2500000);
   dbyte_swap(tmp,2500000);
   dbyte_swap(tmp,2500000);   
   dbyte_swap(tmp,2500000);
   dbyte_swap(tmp,2500000);
   dbyte_swap(tmp,2500000);
   t2=clock();
   dif = difftime(t2, t1);
   printf("old: %lf s \n",dif);
   t1=clock();
   qcd_swap_8(tmp,2500000);
   qcd_swap_8(tmp,2500000);
   qcd_swap_8(tmp,2500000);
   qcd_swap_8(tmp,2500000);
   qcd_swap_8(tmp,2500000);            
   qcd_swap_8(tmp,2500000);               
   t2=clock();
   dif = difftime(t2, t1);
   printf("new: %lf s \n",dif);
   for (i=0; i<2500000; i++)
   {
      tmp[i] -= i+1/(i+1);
      if(tmp[i] != 0) printf("oops!\n");
   }
   free(tmp);
   dif = 0.0;
   dbyte_swap(&dif,1);
   printf("%e\n",dif);
   dif = 0.0;
   qcd_swap_8(&dif,1);
   printf("%e\n",dif);
   printf("%i\n",sizeof(char*));
}
