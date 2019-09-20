// Copyright 2019 SBlincov
#include <ctime>
#include <cmath>
#include <iostream>
#include "omp.h"

#define DEM 100
#define EPS 0.001

float A[DEM][DEM];
float F[DEM];
float Xk[DEM], Zk[DEM];
float Rk[DEM], Sz[DEM], alpha, beta, mf;
float Spr, Spr1, Spz;

float solve(float A[DEM][DEM], float F[DEM]) {
    int i, j, kl = 1;
    int MAX_ITERATIONS = 32767;

    for (mf = 0, i = 0; i < DEM; i++) {
        mf += F[i] * F[i];
    }

    for (i = 0; i < DEM; i++) {
        Xk[i] = 0.2;
    }

    for (i = 0; i < DEM; i++) {
        for (Sz[i]=0, j = 0; j < DEM; j++)
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

    std::cout << "Ответ:"<< std::endl;
    for (i = 0; i < DEM; i++)
        std::cout << "(" << Xk[i] << ")" << std::endl;

    return 0;
}


int main() {
    auto seed = static_cast<unsigned int>(time(0));

    std::cout << "A:" << std::endl;
    for (int i = 0; i < DEM; i++) {
        std::cout << "( ";
        for (int j = 0; j < DEM; j++) {
            A[i][i] = rand_r(&seed) % 100;
            if ((i != j) && (i < j)) {
                A[i][j] = rand_r(&seed) % 100;
            }
            A[j][i] = A[i][j];
            std::cout << A[j][i] << " ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << "b:" << std::endl;
    for (int j = 0; j < DEM; j++) {
        F[j] = rand_r(&seed) % 100;
        std::cout << "( " << F[j] << " )" << std::endl;
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
    solve(A, F);
    time = clock() - time;
    std::cout << "Время выполнения = " << static_cast<double>(time) << "мс" << std::endl;
    return 0;
}
