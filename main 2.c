#include <stdio.h>
#include <math.h>
#include "cmath.h"

#define matrixSize 8
#define baseFirstElement 6

const double E[matrixSize][matrixSize] = {{ 1, 0, 0, 0, 0, 0, 0, 0 },
                                          { 0, 1, 0, 0, 0, 0, 0, 0 },
                                          { 0, 0, 1, 0, 0, 0, 0, 0 },
                                          { 0, 0, 0, 1, 0, 0, 0, 0 },
                                          { 0, 0, 0, 0, 1, 0, 0, 0 },
                                          { 0, 0, 0, 0, 0, 1, 0, 0 },
                                          { 0, 0, 0, 0, 0, 0, 1, 0 },
                                          { 0, 0, 0, 0, 0, 0, 0, 1 }};

void multiply(const double left[matrixSize][matrixSize], const double right[matrixSize][matrixSize], double result[matrixSize][matrixSize])
{
  for(int k = 0; k < matrixSize; k++) {
    for(int i = 0; i < matrixSize; i++) {
      result[k][i] = 0;
      for(int j = 0; j < matrixSize ; j++) {
        result[k][i] +=  left[k][j] * right[j][i];
      }
    }
  }
}

double firstMatrixNorm(double matrix[matrixSize][matrixSize])
{
  double rowSum = 0;
  double result = 0;

  for(int i = 0; i < matrixSize; i++) {
    for(int j = 0; j < matrixSize; j++) {
      matrix[i][j] = fabs(matrix[i][j]);
      rowSum += matrix[i][j];
    }
    if( result < rowSum ) {
      result = rowSum;
    }
  }
  return result;
}

void residualMatrix(const double matrix[matrixSize][matrixSize], double result[matrixSize][matrixSize])
{
  for (int i = 0; i < matrixSize ; i++) {
    for (int j = 0; j < matrixSize; j++) {
      result[i][j] = E[i][j] - matrix[i][j];
    }
  }
}

void inverse(double A[matrixSize][matrixSize], double *cond, double result[matrixSize][matrixSize])
{
  double B[matrixSize];
  double A_copy[matrixSize * matrixSize];
  int ipvt[matrixSize], flag;

  for (int i = 0; i < matrixSize; i++) {
    for (int j = 0; j < matrixSize; j++) {
      A_copy[i*matrixSize + j] = A[i][j];
    }
  }

  decomp(matrixSize, matrixSize, A_copy, cond, ipvt, &flag);

  for (int i = 0; i < matrixSize; i++) {
    for (int j = 0; j < matrixSize; j++) {
      B[j] = E[j][i];
    }
    solve(matrixSize, matrixSize, A_copy, B, ipvt);
    for (int j = 0; j < matrixSize; j++) {
      result[j][i] = B[j];
    }
  }
}

int main() {
  double cond;
  double norm;

  double A[matrixSize][matrixSize] = {{  0,   2,   6,  8, -2,   1,   8, -5},
                                      {  6, -22,  -2, -1,  0,   5,  -6,  4},
                                      { -2,  -3, -16,  0,  0,  -4,   2, -5},
                                      {  1,   1,   4,  9,  1,   0,   0, -6},
                                      {  0,   2,   0,  2, -3,  -5,   7,  5},
                                      {  6,  -2,  -4,  2, -8, -12,   3, -3},
                                      { -6,  -6,   0, -8,  0,   5, -15,  0},
                                      {  0,   7,   6,  0, -5,  -8,  -5, -3}};

  const double p[5] = {1.0, 0.1, 0.01, 0.0001, 0.000001};
  
  double inverseA[matrixSize][matrixSize];
  double AA_1[matrixSize][matrixSize];
  double R[matrixSize][matrixSize];

  for (int i = 0; i < 5; i++) {
    A[0][0] = p[i] + baseFirstElement;

    inverse(A, &cond, inverseA);
    multiply(A, inverseA, AA_1);
    residualMatrix(AA_1, R);
    norm = firstMatrixNorm(R);
    printf("p:    %lf\n", p[i]);
    printf("cond: %lf\n", cond);
    printf("norm: %e\n\n", norm);
  }

  getchar();
  return 0;
}
