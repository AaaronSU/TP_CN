/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "math.h"

void eig_poisson1D(double *eigval, int *la)
{
  for (int i = 1; i < *la + 1; ++i)
  {
    eigval[i - 1] = 2 - (2 * cos((M_PI * i)) / (*la + 1));
  }
}

double eigmax_poisson1D(int *la)
{
  return 2 - (2 * cos((M_PI * 1)) / (*la + 1));
}

double eigmin_poisson1D(int *la)
{
  return 2 - (2 * cos((M_PI * *la)) / (*la + 1));
}

double richardson_alpha_opt(int *la)
{
  return 2 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

// x(k+1) = x(k) + alpha * (b-Ax(k))
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
  // cblas_dgbmv("N", la, lab, kl, ku, alpha_rich, AB, la, X, 1, X, 1);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite)
{
}
