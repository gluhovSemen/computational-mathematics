#include <iostream>
#include "Forsythe.h"
#include <conio.h>
#include <math.h>
#include <iomanip>
#include <fstream>


void f(double x, double y[], double dydx[]) {

    dydx[0] = -4*y[0] + 23*y[1] + exp(-x);
    dydx[1] = 4*y[0] - 48*y[1] + sin(x);
}

void rk3(void(*f)(double x, double y[], double dydx[]),
          double t, double tout, double Zn[],  double h);

int main() {

    double x[2] = {1, 0};

    // Решение с помощью программы rkf45
    int neqn = 2;
    unsigned char work[6*(neqn*sizeof(Float)) + sizeof(struct rkf_inside)];

    rkf ARG;
    ARG.f = f;
    ARG.neqn = neqn;
    ARG.re = 0.0001;
    ARG.ae = 0.0001;
    ARG.work = work;
    ARG.flag = 1;
    ARG.Y = x;
    ARG.t = 0;

    std::ofstream out("out.txt");

    out <<"************* RKF45 *************\n\n";
    out <<"   t      x[0]       x[1]    Flag\n";
    out <<"---------------------------------\n";
    out.setf(std::ios::fixed);

    for(double h=0.1; h<2.1; h+=0.1) {

        ARG.tout = h;
        rkf45(&ARG);
        out << std::setw(5)  << std::setprecision(1) << ARG.t
                  << std::setw(11) << std::setprecision(6) << x[0]
                  << std::setw(11) << std::setprecision(6) << x[1]
                  << std::setw(4)  << std::setprecision(0) << ARG.flag
                  << std::endl;
    }
    out <<"---------------------------------\n\n";


    // Метод Рунге-Кутты 3-й степени точности

    out <<"* RUNGE-KUTTA THIRD-ORDER METHOD *\n\n";
    out.precision(6);
    out.unsetf(std::ios::fixed);

    double t = 0, tout = 2;
    double h;

    h = 0.1;
    out <<"\n h = "<< h << std::endl;
    out <<"   t          x[0]           x[1]   \n";
    out <<"------------------------------------\n";
    for(double tMid = h; tMid < tout + 0.01; tMid += h) {
        x[0] = 1;
        x[1] = 0;
        rk3(f, t, tMid, x, h);
        out << std::setw(5)  << std::setprecision(3) << tMid
            << std::setw(15) << std::setprecision(6) << x[0]
            << std::setw(15) << std::setprecision(6) << x[1]
            << std::endl;
    }

    out.setf(std::ios::fixed);
    h = 0.01;
    out <<"\n h = "<< h << std::endl;
    out <<"   t          x[0]           x[1]   \n";
    out <<"------------------------------------\n";
    int k = 0;
    for(double tMid = h; tMid < tout + 0.01; tMid += 0.01) {
        k++;
        x[0] = 1;
        x[1] = 0;
        rk3(f, t, tMid, x, h);
        if(k==10) {
            out << std::setw(5)  << std::setprecision(1) << tMid
                << std::setw(15) << std::setprecision(6) << x[0]
                << std::setw(15) << std::setprecision(6) << x[1]
                << std::endl;
                k = 0;
        }
    }

    std::cout <<"\nPress any key to continiue...";
    getch();
    return 0;
}

void rk3(void(*f)(double x, double y[], double dydx[]),
          double t, double tout, double Zn[],  double h) {

    double k1[2], k2[2], k3[2];
    double fun[2], work[2];
    double tn;

    for(double t=0; t<=tout; t+=h) {

        tn = t;
        work[0] = Zn[0];
        work[1] = Zn[1];
        f(tn, work, fun);
        k1[0] = h*fun[0];
        k1[1] = h*fun[1];

        tn = t + h/2;
        work[0] = Zn[0] + k1[0]/2;
        work[1] = Zn[1] + k1[1]/2;
        f(tn, work, fun);
        k2[0] = h*fun[0];
        k2[1] = h*fun[1];

        tn = t+h;
        work[0] = Zn[0] - k1[0] + 2*k2[0];
        work[1] = Zn[1] - k1[1] + 2*k2[1];
        f(tn, work, fun);
        k3[0] = h*fun[0];
        k3[1] = h*fun[1];

        Zn[0] = Zn[0] + (k1[0] + 4*k2[0] + k3[0]) / 6;
        Zn[1] = Zn[1] + (k1[1] + 4*k2[1] + k3[1]) / 6;
    }
}

