#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "eigen_algorithms.h"

void tridiagonalize_symmetric(int n, double *A, double *diag, double *subdiag) {
    int i, j, k;
    double x, y, r, c, s;
    double tmp1, tmp2;
    
    for (j = 0; j < n - 2; j++) {
        for (i = n - 1; i > j + 1; i--) {
            x = A[(i-1) * n + j];
            y = A[i * n + j];
            
            if (fabs(y) < 1e-15) continue;
            
            r = sqrt(x * x + y * y);
            if (r < 1e-15) continue;
            
            c = x / r;
            s = -y / r;
            
            for (k = 0; k < n; k++) {
                tmp1 = A[(i-1) * n + k];
                tmp2 = A[i * n + k];
                A[(i-1) * n + k] = c * tmp1 - s * tmp2;
                A[i * n + k] = s * tmp1 + c * tmp2;
            }
            
            for (k = 0; k < n; k++) {
                tmp1 = A[k * n + (i-1)];
                tmp2 = A[k * n + i];
                A[k * n + (i-1)] = c * tmp1 - s * tmp2;
                A[k * n + i] = s * tmp1 + c * tmp2;
            }
        }
    }
    
    for (i = 0; i < n; i++) {
        diag[i] = A[i * n + i];
        if (i < n - 1) {
            subdiag[i] = A[i * n + (i+1)];
        }
    }
}

int sturm_count(int n, double *diag, double *subdiag, double lambda) {
    double p0 = 1.0;
    double p1 = diag[0] - lambda;
    int count = (p1 < 0) ? 1 : 0;
    int i;
    
    for (i = 1; i < n; i++) {
        double p2 = (diag[i] - lambda) * p1 - subdiag[i-1] * subdiag[i-1] * p0;
        
        if ((p1 < 0 && p2 >= 0) || (p1 >= 0 && p2 < 0)) {
            count++;
        }
        
        p0 = p1;
        p1 = p2;
        
        if (fabs(p1) > 1e100) {
            p0 /= 1e100;
            p1 /= 1e100;
        }
    }
    
    return count;
}

double bisection_kth_eigenvalue(int n, double *diag, double *subdiag, int k, double eps, int *iter_count) {
    double left = 1e100;
    double right = -1e100;
    double mid, center, radius;
    int target, count, i;
    
    for (i = 0; i < n; i++) {
        radius = 0.0;
        if (i > 0) radius += fabs(subdiag[i-1]);
        if (i < n - 1) radius += fabs(subdiag[i]);
        
        center = diag[i];
        if (center - radius < left) left = center - radius;
        if (center + radius > right) right = center + radius;
    }
    
    left -= 0.1 * (right - left);
    right += 0.1 * (right - left);
    
    target = n - k + 1;
    
    *iter_count = 0;
    
    while (right - left > eps && *iter_count < 10000) {
        mid = (left + right) / 2.0;
        count = sturm_count(n, diag, subdiag, mid);
        
        if (count >= target) {
            right = mid;
        } else {
            left = mid;
        }
        
        (*iter_count)++;
    }
    
    return (left + right) / 2.0;
}