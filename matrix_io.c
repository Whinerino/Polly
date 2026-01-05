#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_io.h"

double formula(int k, int n, int i, int j) {
    int ii = i + 1;
    int jj = j + 1;
    
    switch (k) {
        case 1:
            return n - ((ii > jj) ? ii : jj) + 1.0;
        case 2:
            if (ii == jj) return 2.0;
            if (abs(ii - jj) == 1) return -1.0;
            return 0.0;
        case 3:
            if (ii == jj && ii < n) return 1.0;
            if (jj == n) return (double)ii;
            if (ii == n) return (double)jj;
            return 0.0;
        case 4:
            return 1.0 / (ii + jj - 1);
        default:
            return 0.0;
    }
}

void read_matrix(int n, double *A, int k, char *filename) {
    int i, j;
    FILE *file;
    double avg;
    
    if (k == 0) {
        file = fopen(filename, "r");
        if (!file) {
            printf("Error: Ñannot open file %s\n", filename);
            exit(1);
        }
        
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (fscanf(file, "%lf", &A[i * n + j]) != 1) {
                    printf("Error: Invalid data format in file\n");
                    fclose(file);
                   exit(1);
                }
            }
        }
        fclose(file);
    } else {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                A[i * n + j] = formula(k, n, i, j);
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            avg = (A[i * n + j] + A[j * n + i]) / 2.0;
            A[i * n + j] = avg;
            A[j * n + i] = avg;
        }
    }
}

void print_matrix(int n, int m, double *A) {
    int limit = (m < n) ? m : n;
    int i, j;
    
    for (i = 0; i < limit; i++) {
        for (j = 0; j < limit; j++) {
            printf(" %10.3e", A[i * n + j]);
        }
        printf("\n");
    }
	printf("\n");
}