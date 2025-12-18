#include "matrix_op.h"
#include <stdio.h>

static double det2(double a, double b, double c, double d) {
    return a * d - b * c;
}

void matrix_add(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            R[i][j] = A[i][j] + B[i][j];
}

void matrix_sub(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            R[i][j] = A[i][j] - B[i][j];
}

void matrix_hadamard(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            R[i][j] = A[i][j] * B[i][j];
}

void matrix_mul(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            double sum = 0.0;
            for (int k = 0; k < SIZE; k++) {
                sum += A[i][k] * B[k][j];
            }
            R[i][j] = sum;
        }
    }
}

void matrix_transpose(const double A[SIZE][SIZE], double R[SIZE][SIZE]) {
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            R[j][i] = A[i][j];
}

double matrix_det(const double A[SIZE][SIZE]) {
    // |a b c|
    // |d e f| = a(ei - fh) - b(di - fg) + c(dh - eg)
    // |g h i|
    double a = A[0][0], b = A[0][1], c = A[0][2];
    double d = A[1][0], e = A[1][1], f = A[1][2];
    double g = A[2][0], h = A[2][1], i = A[2][2];

    return a * det2(e, f, h, i)
         - b * det2(d, f, g, i)
         + c * det2(d, e, g, h);
}

void matrix_adjoint(const double A[SIZE][SIZE], double R[SIZE][SIZE]) {
    // adj(A) = cofactor(A)^T
    // Cofactor C_ij = (-1)^(i+j) * det(M_ij)
    // Here build cofactors then transpose into R.

    double C[SIZE][SIZE];

    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            // build 2x2 minor matrix excluding row i and col j
            double m[2][2];
            int r = 0;
            for (int rr = 0; rr < SIZE; rr++) {
                if (rr == i) continue;
                int c = 0;
                for (int cc = 0; cc < SIZE; cc++) {
                    if (cc == j) continue;
                    m[r][c] = A[rr][cc];
                    c++;
                }
                r++;
            }
            double minor_det = det2(m[0][0], m[0][1], m[1][0], m[1][1]);
            double sign = ((i + j) % 2 == 0) ? 1.0 : -1.0;
            C[i][j] = sign * minor_det;
        }
    }

    // transpose cofactors -> adjoint
    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            R[i][j] = C[j][i];
}

int matrix_inverse(const double A[SIZE][SIZE], double R[SIZE][SIZE]) {
    double detA = matrix_det(A);
    if (detA == 0.0) return 0; // singular

    double adj[SIZE][SIZE];
    matrix_adjoint(A, adj);

    for (int i = 0; i < SIZE; i++)
        for (int j = 0; j < SIZE; j++)
            R[i][j] = adj[i][j] / detA;

    return 1;
}

void matrix_print(const char *name, const double A[SIZE][SIZE]) {
    if (name) printf("%s =\n", name);
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            printf("%10.4f ", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}
