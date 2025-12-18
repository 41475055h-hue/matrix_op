#include <stdio.h>
#include "matrix_op.h"

int main(void) {
    double A[SIZE][SIZE] = {
        {1, 2, 3},
        {0, 1, 4},
        {5, 6, 0}
    };

    double B[SIZE][SIZE] = {
        { -1,  0,  2},
        {  3,  1,  0},
        {  4, -2,  1}
    };

    double R[SIZE][SIZE];

    matrix_print("A", A);
    matrix_print("B", B);

    matrix_add(A, B, R);
    matrix_print("A + B", R);

    matrix_sub(A, B, R);
    matrix_print("A - B", R);

    matrix_hadamard(A, B, R);
    matrix_print("A o B (element-wise)", R);

    matrix_mul(A, B, R);
    matrix_print("A * B", R);

    matrix_transpose(A, R);
    matrix_print("A^T", R);

    double detA = matrix_det(A);
    printf("det(A) = %.4f\n\n", detA);

    matrix_adjoint(A, R);
    matrix_print("adj(A)", R);

    if (matrix_inverse(A, R)) {
        matrix_print("A^-1", R);
    } else {
        printf("A is singular, no inverse.\n");
    }

    return 0;
}
