#include <stdio.h>
#include <iostream>
#include <cmath>
#include <math.h>
using namespace std;
double x0 = 0.525; // точка разрыва

void progonka(double* x, double* y, double* k, double* q, double* fi, int N, int m) {
    double* alfa; // прогоночные коэффиценты
    double* beta;
    double* a; // диагональные элементы матрицы
    double* b;
    double* c;

    double h = 1.0 / N; // шаг

    //cout << "h =" << h << '\n';

    for (int i = 0; i < N + 1; i++) {
        x[i] = 0.0;
        y[i] = 0.0;
        k[i] = 0.0;
        q[i] = 0.0;
        fi[i] = 0.0;
    }
    for (int i = 0; i < N + 1; i++) {
        x[i] = i * h;
    }

    for (int i = 0; i < N + 1; i++) {
        if (m == 1) { // немодельная задача
            if (x[i] < x0) {
                k[i] = exp(-x[i] * x[i]);
                q[i] = x[i] * x[i];
                fi[i] = sin(x[i]);
            }
            else {
                k[i] = x[i];
                q[i] = x[i] * x[i];
                fi[i] = sin(x[i]);
            }
        } else { // модельная задача
            if (x[i] < x0) {
                k[i] = exp(-x0 * x0);
                q[i] = x0 * x0;
                fi[i] = sin(x0);
            } else {
                k[i] = x0;
                q[i] = x0 * x0;
                fi[i] = sin(x0);
            }
        }
    }

    alfa = new double[N + 2];
    beta = new double[N + 2];
    a = new double[N + 2];
    b = new double[N + 2];
    c = new double[N + 2];

    

    for (int i = 1; i < N; i++) {
        a[i] = k[i] / (h * h);
        //cout << "a[" << i << "]=" << a[i] << " " << '\n';
        b[i] = k[i + 1] / (h * h);
        //cout << "b[" << i << "]=" << b[i] << " " << '\n';
        c[i] = k[i] / (h * h) + k[i + 1] / (h * h) + q[i];
        //cout << "c[" << i << "]=" << c[i] << " " << '\n';
    }
    b[0] = 0;
    a[N] = 0;
    c[0] = c[N] = 1;
    fi[0] = 0; fi[N] = 1;
    alfa[1] = 0.0;
    beta[1] = 0.0;

    // прямой ход прогонки
    for (int i = 1; i < N; i++) {
        alfa[i + 1] = b[i] / (c[i] - a[i] * alfa[i]);
        beta[i + 1] = (fi[i] + a[i] * beta[i]) / (c[i] - alfa[i] * a[i]);
    }
    // обратный ход прогонки
    y[N] = 1;
    for (int i = N - 1; i >= 0; i--) {
        y[i] = alfa[i + 1] * y[i + 1] + beta[i + 1];
    }
    delete[]a;
    delete[]b;
    delete[]c;
    delete[]alfa;
    delete[]beta;
}
void real_decision(double* u, int N) {
    double h = 1.0 / N;
    double* x;
    x = new double[N + 1];
    for (int i = 1; i < N + 1; i++) {
        x[i] = i * h;
    }
    double k1 = exp(-x0 * x0);
    double q1 = x0 * x0;
    double fi1 = sin(x0);
    double d1 = fi1 / q1; // частное решение

    double k2 = x0;
    double q2 = x0 * x0;
    double fi2 = sin(x0);
    double d2 = fi2 / q2;

    double lambda1 = sqrt(q1 / k1);
    double lambda2 = sqrt(q2 / k2);

    double A1 = exp(-lambda1 * x0) - exp(lambda1 * x0);
    double B1 = exp(-2 * lambda2 + lambda2 * x0) - exp(-lambda2 * x0);
    double D1 = d1 * (exp(lambda1 * x0) - 1) + exp(lambda2 * x0 - lambda2) * (1 - d2) + d2;

    double A2 = -k1 * (lambda1 * exp(lambda1 * x0) + lambda1 * exp(-lambda1 * x0));
    double B2 = k2 * (lambda2 * exp(-2 * lambda2 + lambda2 * x0) + lambda2 * exp(-lambda2 * x0));
    double D2 = k1 * d1 * lambda1 * exp(lambda1 * x0) + k2 * lambda2 * exp(-lambda2 + lambda2 * x0) * (1 - d2);

    double znam = A1 * B2 - B1 * A2;

    double const2 = (D1 * B2 - B1 * D2) / znam;
    double const4 = (A1 * D2 - A2 * D1) / znam;
    double const1 = -d1 - const2;
    double const3 = exp(-lambda2) * (1 - d2 - const4 * exp(-lambda2));
    

    u[0] = 0;
    u[N] = 1;
    for (int i = 1; i < N; i++) {
        if (x[i] < x0)
            u[i] = const1 * exp(lambda1 * x[i]) + const2 * exp(-lambda1 * x[i]) + d1;
        else
            u[i] = const3 * exp(lambda2 * x[i]) + const4 * exp(-lambda2 * x[i]) + d2;
    }
    delete[] x;
}
double max_delta(double* u1, double* u2, int N) {
    double Max = 0.0;
    for (int i = 0; i < N + 1; i++) {
        if (Max < abs(u1[i] - u2[i])) {
            Max = abs(u1[i] - u2[i]);
        }
    }

    return Max;
}

int main() {
    int N, m;
    double delta, maximum;
    double* q, * k, * y, * u, * x, * fi;
    cout << "Start program! Enter N, please:" << endl;
    cin >> N;
    double h = 1.0 / N; // шаг
    cout << "Model or no ? Enter m,please: " << endl;
    cin >> m;
    if (m == 0) {
        cout << " model " << endl;
    }
    else {
        cout << "no model" << endl;
    }
    while (N < 10000) {
        if (m == 0) {
            u = new double[N + 1];
            real_decision(u, N);
            x = new double[N + 1];
            y = new double[N + 1];
            q = new double[N + 1];
            k = new double[N + 1];
            fi = new double[N + 1];
            progonka(x, y, k, q, fi, N,m);
            for (int i = 0; i < N + 1; i++) {
                cout <<"x="<< x[i]<< "  u[" << x[i] << "] =" << u[i] << "  y["<< x[i]<<"] =" << y[i] << "  delta = " << abs(y[i] - u[i]) << '\n';
            }
            delta = max_delta(u, y, N);
            if (max_delta(u, y, N) < 0.01) {
                cout << "max delta = " << max_delta(u, y, N) << " < 0.01 " << '\n';
                break;
            }
            else {
                cout << "max delta = " << max_delta(u, y, N) << " > 0.01 " << '\n';
                N = 2 * N;
                cout << "New N = " << N << '\n';
            }
            delete[] x;
            delete[] k;
            delete[] q;
            delete[] fi;
            delete[] y;
            delete[] u;

        }  else {
            maximum = 0;
            x = new double[N + 1];
            y = new double[N + 1];
            q = new double[N + 1];
            k = new double[N + 1];
            fi = new double[N + 1];

            progonka(x, y, k, q, fi, N,m);
            delete[] x;
            delete[] k;
            delete[] q;
            delete[] fi;

            //N = N * 2;
            x = new double[2*N + 1];
            q = new double[2*N + 1];
            k = new double[2*N + 1];
            fi = new double[2*N + 1];
            u = new double[2*N + 1];
            progonka(x, u, k, q, fi, 2*N,m);
            for (int i = 0; i < N + 1; i++) {
                cout <<" x = " << x[2*i] <<" y_N[" << x[2 * i] << "] = " << y[i] << " y_2N[" << x[2 * i] << "] = " << u[2 * i] << " delta = " << abs(u[2 * i] - y[i]) << '\n';
                if (maximum < abs(u[2 * i] - y[i]))
                    maximum = abs(u[2 * i] -y[i]);
            }
            if (maximum < 0.01) {
                cout << "max delta = " << maximum << " < 0.01 " << '\n';
                break;
            }
            else {
                cout << "max delta = " << maximum << " > 0.01 " << '\n';
                N = 2*N;
                cout << "New N = " << N << '\n';
            }

            delete[] x;
            delete[] k;
            delete[] q;
            delete[] fi;
            delete[] y;
            delete[] u;
        } 
    }
    return 0;
}