/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, int *kv)
{
  int rl = *la;
  int cl = *lab;
  int col0 = *kv;
  int i = 1;
  while (i < rl * cl - 1)
  {
    if (i % cl == 0)
    {
      for (int j = 0; j < col0; ++j)
      {
        AB[i] = 0;
        i++;
      }
    }
    if (i == col0)
    {
      AB[i] = 0;
    }
    else if ((i % cl == col0) || (i % cl == col0 + 2))
    {
      AB[i] = -1;
    }
    else if (i % cl == col0 + 1)
    {
      AB[i] = 2;
    }
    i++;
  }
  AB[rl * cl - 1] = 0;
}

void set_GB_operator_colMajor_poisson1D_Id(double *AB, int *lab, int *la, int *kv)
{
  int rl = *la;
  int cl = *lab;
  int col0 = *kv;
  int i = 1;
  while (i < rl * cl - 1)
  {
    if (i % cl == col0 + 1)
    {
      AB[i] = 1;
    }
    else
    {
      AB[i] = 0;
    }
    i++;
  }
  AB[rl * cl - 1] = 0;
}

void set_dense_RHS_DBC_1D(double *RHS, int *la, double *BC0, double *BC1)
{
  RHS[0] = *BC0;
  for (int i = 1; i < *la - 1; ++i)
  {
    RHS[i] = 0;
  }
  RHS[*la - 1] = *BC1;
}

void set_analytical_solution_DBC_1D(double *EX_SOL, double *X, int *la, double *BC0, double *BC1)
{
  const double DT = *BC1 - *BC0;

  for (int i = 0; i < (*la); i++)
  {
    EX_SOL[i] = *BC0 + X[i] * DT; // Equation analytique
  }
}

void set_grid_points_1D(double *x, int *la)
{
  int jj;
  double h;
  h = 1.0 / (1.0 * ((*la) + 1));
  for (jj = 0; jj < (*la); jj++)
  {
    x[jj] = (jj + 1) * h;
  }
}

double relative_forward_error(double *x, double *y, int *la)
{
  cblas_daxpy(*la, -1, x, 1, y, 1);
  double sum_bas = cblas_dnrm2(*la, x, 1);
  double sum_haut = cblas_dnrm2(*la, y, 1);

  // En théorie, nous devons trouver 10e-16, en pratique, ceci n'est pas le cas. A voir.
  // Déjà résolu, l'erreur venant de la fonction void set_grid_points_1D(double *x, int *la).
  // La fonction donne des résultats précis mais de l'ordre de 10^-8
  return sum_haut / sum_bas;
}

int indexABCol(int i, int j, int *lab)
{
  return j * (*lab) + i;
}

int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info)
{
  for (int ii = *lab - 1; ii < *lab * *la; ii++)
  {
    if (!(ii % *lab))
    {
      if (AB[ii - 2])
      {
        AB[ii - 1] = AB[ii + 1] / AB[ii - 2];
      }
      AB[ii + 2] = AB[ii + 2] - AB[ii - 1] * AB[ii + 1];
    }
  }

  for (int ip = 0; ip < *la; ip++)
  {
    ipiv[ip] = ip + 1;
  }
  *info = 0;

  return 0;
}