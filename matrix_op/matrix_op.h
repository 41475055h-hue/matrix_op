#ifndef MATRIX_OP_H
#define MATRIX_OP_H

#define SIZE 3

// Basic
void matrix_add(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]);
void matrix_sub(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]);
void matrix_hadamard(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]); // element-wise

// Linear
void matrix_mul(const double A[SIZE][SIZE], const double B[SIZE][SIZE], double R[SIZE][SIZE]);
void matrix_transpose(const double A[SIZE][SIZE], double R[SIZE][SIZE]);

// Advanced
double matrix_det(const double A[SIZE][SIZE]);
void matrix_adjoint(const double A[SIZE][SIZE], double R[SIZE][SIZE]);

// Inverse
// return 1 if invertible, 0 if singular
int matrix_inverse(const double A[SIZE][SIZE], double R[SIZE][SIZE]);

// Utilities (optional but convenient)
void matrix_print(const char *name, const double A[SIZE][SIZE]);

#endif
