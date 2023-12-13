/******************************************/
/* tp_env.c                               */
/* This file contains a main function to  */
/* test the environment of compilation    */
/******************************************/
#include "tp_env.h"
int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  printf("--------- Test environment of execution for Practical exercises of Numerical Algorithmics ---------\n\n");
  printf("The exponantial value is e = %f \n",M_E);
  printf("The maximum single precision value from values.h is maxfloat = %e \n",MAXFLOAT);
  printf("The maximum single precision value from float.h is flt_max = %e \n",FLT_MAX);
  printf("The maximum double precision value from float.h is dbl_max = %e \n",DBL_MAX);
  printf("The epsilon in single precision value from float.h is flt_epsilon = %e \n",FLT_EPSILON);
  printf("The epsilon in double precision value from float.h is dbl_epsilon = %e \n",DBL_EPSILON);

  printf("\n\n Test of ATLAS (BLAS/LAPACK) environment \n");

  double x[5], y[5];
  int ii;
  for (ii=0;ii<5;ii++){
    x[ii]=ii+1; y[ii]=ii+6;
    printf("x[%d] = %lf, y[%d] = %lf\n",ii,x[ii],ii,y[ii]);
  }

  printf("\nTest DCOPY y <- x \n");
  cblas_dcopy(5,x,1,y,1);
  for (ii=0;ii<5;ii++){
    printf("y[%d] = %lf\n",ii,y[ii]);
  }

  printf("\n\n--------- End -----------\n");
}
