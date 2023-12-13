/**********************************************/
/* lib_poisson1D.h                            */
/* Header for Numerical library developed to  */ 
/* solve 1D Poisson problem (Heat equation)   */
/**********************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int* lab, int *la, int *kv);
void set_GB_operator_colMajor_poisson1D_Id(double* AB, int* lab, int *la, int *kv);
void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1);
void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1);
void set_grid_points_1D(double* x, int* la);
void write_GB2AIJ_operator_poisson1D(double* AB, int *la, char* filename);
void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int *la, char* filename);
void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename);
void write_vec(double* vec, int* la, char* filename);
void write_xy(double* vec, double* x, int* la, char* filename);
int indexABCol(int i, int j, int *lab);
int dgbtrftridiag(int *la, int *n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info);
