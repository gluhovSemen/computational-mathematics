#include <stdio.h>
#include <math.h>
#include "cmath.h"


const int lowerBound = 0;
#define upperBound 1
#define lagrangDegree 6
#define pointsNumber 7
#define pointsApprox 6
#define minLimit 0
#define maxLimit 3
#define step 0.5  

double x;

double function(double z)
{
  return 1/(exp(1.9 * z * z) + x);
}

double lagrang(double pointX, double* x, double* y, int n)
{
  double lagrSum = 0;

  for (int i = 0; i < n; i++) {
    double P = 1;
    for (int j = 0; j < n; j++) {
      if (i != j) {
        P *= (pointX - x[j]) / (x[i] - x[j]);
      }
    }
    lagrSum += P * y[i];
  }
  return lagrSum;
}

int main() { 
  double abserr = 0.00000001;
  double relerr = 0.00000001;
  
  double res;

  double *errestR = NULL;
  int nofunR = 0;
  double *posnR = NULL;
  int *flag = NULL;
  
  double x0[pointsNumber];
  double y0[pointsNumber];

 
  for (int count = 0; count < pointsNumber; count++) {
    x = 0.5 * count;
    x0[count] = x;
    quanc8(&function, lowerBound, upperBound, abserr, relerr, &res, &errestR, &nofunR, &posnR, &flag);
    y0[count] = res;
  }

  double l[pointsNumber], s[pointsNumber];

  double B[pointsNumber], C[pointsNumber], D[pointsNumber];
  int flagS = 0;
  int last = 0;
  int e1 = 1;
  int e2 = 0;
  double s1 = 0;
  double s2 = 0;

  double xk[pointsApprox];
  double yk[pointsApprox];

  for (int i = 0; i < pointsApprox; i++) {
    xk[i] = 0.25 + 0.5 * i;
    x = xk[i];
    quanc8(&function, lowerBound, upperBound, abserr, relerr, &res, &errestR, &nofunR, &posnR, &flag);
    yk[i] = res;
  }

  spline(pointsApprox, e1, e2, s1, s2, xk, yk, B, C, D, &flagS);
  for (int k = 0; k < pointsApprox; k++) {
    l[k] = lagrang(xk[k], x0, y0, lagrangDegree);
    s[k] = seval(pointsApprox, xk[k], x0, y0, B, C, D, &last);
  }

  for (int k = 0; k < pointsApprox; k++) {
    printf("%lf                \n %lf        \n%lf  \n          %lf   \n", xk[k],yk[k], l[k], s[k]);
  }

  getchar();
  return 0;
}
