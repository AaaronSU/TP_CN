/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void write_GB_operator_rowMajor_poisson1D(double *AB, int *lab, int *la, char *filename)
{
  FILE *file;
  int ii, jj;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL)
  {
    for (ii = 0; ii < (*lab); ii++)
    {
      for (jj = 0; jj < (*la); jj++)
      {
        fprintf(file, "%lf\t", AB[ii * (*la) + jj]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double *AB, int *lab, int *la, char *filename)
{
  FILE *file;
  int ii, jj;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL)
  {
    for (ii = 0; ii < (*la); ii++)
    {
      for (jj = 0; jj < (*lab); jj++)
      {
        fprintf(file, "%lf\t", AB[ii * (*lab) + jj]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double *AB, int *la, char *filename)
{
  FILE *file;
  int jj;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL)
  {
    for (jj = 1; jj < (*la); jj++)
    {
      fprintf(file, "%d\t%d\t%lf\n", jj, jj + 1, AB[(*la) + jj]);
    }
    for (jj = 0; jj < (*la); jj++)
    {
      fprintf(file, "%d\t%d\t%lf\n", jj + 1, jj + 1, AB[2 * (*la) + jj]);
    }
    for (jj = 0; jj < (*la) - 1; jj++)
    {
      fprintf(file, "%d\t%d\t%lf\n", jj + 2, jj + 1, AB[3 * (*la) + jj]);
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}

void write_vec(double *vec, int *la, char *filename)
{
  int jj;
  FILE *file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL)
  {
    for (jj = 0; jj < (*la); jj++)
    {
      fprintf(file, "%lf\n", vec[jj]);
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}

void write_xy(double *vec, double *x, int *la, char *filename)
{
  int jj;
  FILE *file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL)
  {
    for (jj = 0; jj < (*la); jj++)
    {
      fprintf(file, "%lf\t%lf\n", x[jj], vec[jj]);
    }
    fclose(file);
  }
  else
  {
    perror(filename);
  }
}
