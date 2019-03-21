//
// Created by Blintsov Sergey on 10/03/2019.
//
#include <ctime>
#include <cmath>
#include <iostream>
using namespace std;

#define DEM 3
#define EPS 0.001

float A[DEM][DEM];
float F[DEM];
float Xk[DEM], Zk[DEM];
float Rk[DEM], Sz[DEM], alpha, beta, mf;
float Spr, Spr1, Spz;

float solve(float A[DEM][DEM], float F[DEM]) {
    int i,j,kl=1;
    int MAX_ITERATIONS = 32767;

    for (mf = 0,i = 0; i < DEM; i++) {
        mf += F[i] * F[i];
    }

    for (i = 0; i < DEM; i++) {
        Xk[i] = 0.2;
    }

    for (i = 0; i < DEM; i++) {
        for (Sz[i]=0,j = 0; j < DEM; j++)
            Sz[i] += A[i][j] * Xk[j];
        Rk[i] = F[i] - Sz[i];
        Zk[i] = Rk[i];
    }

    int iterationCounter = 0;

    do {
        iterationCounter++;
        Spz = 0;
        Spr = 0;

        for (i = 0; i < DEM; i++) {
            for (Sz[i] = 0, j = 0; j < DEM; j++) {
                Sz[i] += A[i][j] * Zk[j];
            }
            Spz += Sz[i] * Zk[i];
            Spr += Rk[i] * Rk[i];
        }
        alpha = Spr / Spz;

        Spr1 = 0;
        for (i = 0; i < DEM; i++) {
            Xk[i] += alpha * Zk[i];
            Rk[i] -= alpha * Sz[i];
            Spr1 += Rk[i] * Rk[i];
        }
        kl++;

        beta = Spr1 / Spr;

        for (i = 0; i < DEM; i++)
            Zk[i] = Rk[i] + beta * Zk[i];

    } while (Spr1 / mf > EPS * EPS && iterationCounter < MAX_ITERATIONS);

    cout << "Ответ:"<< endl;
    for (i = 0; i < DEM; i++)
        cout << "(" << Xk[i] << ")" << endl;

    return 0;
}


int main() {
    srand(time(0));

    cout << "A:" << endl;
    for (int i = 0; i < DEM; i++) {
        cout << "( ";
        for (int j = 0; j < DEM; j++) {
            A[i][i] = rand() % 2;
            if ((i != j) && (i < j)) {
                A[i][j] = rand() % 2 - rand() % 2;
            }
            A[j][i] = A[i][j];
            cout << A[j][i] << " ";
        }
        cout << ")" << endl;
    }
    cout << "b:" << endl;
    for (int j = 0; j < DEM; j++){
        F[j] = rand() % 2 - rand() % 2;
        cout << "( " << F[j] << " )" << endl;
    }
// Задаем вручную:
//    A[0][0] = 3;
//    A[0][1] = 4;
//    A[0][2] = 0;
//    A[1][0] = 4;
//    A[1][1] = -3;
//    A[1][2] = 0;
//    A[2][0] = 0;
//    A[2][1] = 0;
//    A[2][2] = 5;
//
//    F[0] = 1;
//    F[1] = 5;
//    F[2] = 9;

    clock_t time = clock();
    solve(A,F);
    time = clock() - time;
    cout << "Время выполнения = " << static_cast<double>(time) << "мс" << endl;
    return 0;
}