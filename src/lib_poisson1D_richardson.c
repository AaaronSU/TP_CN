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
  cblas_dcopy(*la, RHS, 1.0, resvec, 1.0);
  cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, resvec, 1);

  double norme2RHS = 1 / cblas_dnrm2(*la, RHS, 1);
  double residu = cblas_dnrm2(*la, resvec, 1) * norme2RHS; // erreur residuelle

  // boucler si le résidu est plus grand que la tolérance et le nombre d'itération est plus petit que le maximum des itérations.
  while (residu > (*tol) && *nbite < *maxit)
  {
    cblas_daxpy(*la, *alpha_rich, resvec, 1, X, 1); // X = X + alpha * (b - Ax)

    cblas_dcopy(*la, RHS, 1.0, resvec, 1.0);
    cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, resvec, 1);

    // Nouvelle erreur résiduelle
    residu = cblas_dnrm2(*la, resvec, 1) * norme2RHS;
    resvec[*nbite] = residu;

    // on incrémente le nombre de boucle
    *nbite += 1;
  }

  printf("\nL'erreur résiduelle obtenu à la fin est de %f\n", residu);
  printf("Le nombre total d'itération pour richardson alpha est de %d\n", *nbite);
}

void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
  // Selon le calcul, la valeur située sur la diagonale est de 0.5.
  for (int ii = 0; ii < *lab * *la; ii++)
  {
    if (ii % *lab == *kl + *kv - 1)
    {
      MB[ii] = AB[ii];
    }
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la, int *ku, int *kl, int *kv)
{
  // Selon le calcul, la valeur située sur la diagonale est de 0.5, tandis que les valeurs situées en dessous de la diagonale sont de 0.25.
  for (int ii = 0; ii < *lab * *la; ii++)
  {
    if (ii % *lab == *kl + *kv - 1 || ii % *lab == *kl + *kv)
    {
      MB[ii] = AB[ii];
    }
  }
}

void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la, int *ku, int *kl, double *tol, int *maxit, double *resvec, int *nbite, int *NRHS, int *ipiv, int *info)
{
  int ku_minus = *ku - 1;
  double norme2RHS = 1 / cblas_dnrm2(*la, RHS, 1);
  double residu = cblas_dnrm2(*la, resvec, 1) * norme2RHS; // erreur residuelle

  // La factorisation LU pour la matrice bande MB
  dgbtrf_(la, la, kl, &ku_minus, MB, lab, ipiv, info);
  // boucler si le résidu est plus grand que la tolérance et le nombre d'itération est plus petit que le maximum des itérations.
  while (*nbite < *maxit)
  {
    cblas_dcopy(*la, RHS, 1.0, resvec, 1.0);
    // b - A * x**k
    cblas_dgbmv(CblasColMajor, CblasConjNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X, 1, 1.0, resvec, 1);

    // Nouvelle erreur résiduelle
    residu = cblas_dnrm2(*la, resvec, 1) * norme2RHS;
    resvec[*nbite] = residu;

    // on incrémente le nombre de boucle
    *nbite += 1;

    dgbtrs_("N", la, kl, &ku_minus, NRHS, MB, lab, ipiv, resvec, la, info, (unsigned long)&la);
    cblas_daxpy(*la, 1, resvec, 1, X, 1);
    if (residu < (*tol))
    {
      break;
    }
  }
  printf("\nL'erreur résiduelle obtenu à la fin est de %f\n", residu);
  printf("Le nombre total d'itération pour richardson MB est de %d", *nbite);
}
