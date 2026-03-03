#ifndef EIGEN_ALGORITHMS_H
#define EIGEN_ALGORITHMS_H

void tridiagonalize_symmetric(int n, double *A, double *diag, double *subdiag);
int sturm_count(int n, double *diag, double *subdiag, double lambda);
double bisection_kth_eigenvalue(int n, double *diag, double *subdiag, int k, double eps, int *iter_count);

#endif