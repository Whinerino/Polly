#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "matrix_io.h"
#include "eigen_algorithms.h"

int main(int argc, char *argv[]) 
{
    int n, m, k_formula, k_eigen, bisect_iterations;
    double eps, t1, t2, lambda_k;
    double *A, *diag, *subdiag;
    char *filename;
    clock_t start_time, tri_time, bisect_time;
    
    if (argc < 5) {
        fprintf(stderr, "Usage: %s n m eps k [filename]\n", argv[0]);
        return 1;
    }

    n = atoi(argv[1]);
    m = atoi(argv[2]);
    eps = atof(argv[3]);
    k_formula = atoi(argv[4]);
    filename = NULL;
    
    if (k_formula == 0) {
        if (argc < 6) {
            fprintf(stderr, "Error: filename required when k=0\n");
            return 1;
        }
        filename = argv[5];
    }

    if (n <= 0 || m <= 0 || eps <= 0 || k_formula < 0 || k_formula > 4) {
        fprintf(stderr, "Error: invalid parameters\n");
        return 1;
    }

    A = (double*)malloc((size_t)(n * n) * sizeof(double));
    diag = (double*)malloc((size_t)n * sizeof(double));
    if (n > 1) {
        subdiag = (double*)malloc((size_t)(n-1) * sizeof(double));
    } else {
        subdiag = NULL;
    }
    
    if (A == NULL || diag == NULL || (n > 1 && subdiag == NULL)) {
        fprintf(stderr, "Memory allocation failed\n");
        free(A);
        free(diag);
        free(subdiag);
        return 1;
    }

    read_matrix(n, A, k_formula, filename);
    
    print_matrix(n, m, A);

    start_time = clock();
    
    tridiagonalize_symmetric(n, A, diag, subdiag);
    
    tri_time = clock();
    t1 = (double)(tri_time - start_time) / CLOCKS_PER_SEC;

    k_eigen = (n / 2) + 1;
    if (k_eigen > n) k_eigen = n;
    if (k_eigen < 1) k_eigen = 1;
    
    bisect_iterations = 0;
    lambda_k = bisection_kth_eigenvalue(n, diag, subdiag, k_eigen, eps, &bisect_iterations);
    
    bisect_time = clock();
    t2 = (double)(bisect_time - tri_time) / CLOCKS_PER_SEC;

    printf("%s : Iterations = %d Iterations1 = %d Elapsed1 = %.2f Elapsed2 = %.2f\n",
           argv[0], bisect_iterations, 
           (n > 0) ? (bisect_iterations + n - 1) / n : 0, t1, t2);
    
    printf("\nResult: %d-th eigenvalue = %.10e\n", k_eigen, lambda_k);

    free(A); 
    free(diag); 
    free(subdiag); 
    
    return 0;
}