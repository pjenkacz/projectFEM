//
// Created by Michał Pieniek on 07/02/2025.
//

#ifndef HEADER_H
#define HEADER_H
#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <sstream>
#include <iomanip>
#include <string>
#include <algorithm>

// npc=4 2punktowy schemat, npc=9 3punktowy schemat, npc=16 4punktowy schemat
#define npc 16
//npcbc=2 2punktowy schemat, npcbc=3 3punktowy schemat, npcbc=4 4punktowy schemat 
#define NPCBC 4
//aby wybrac rodzaj testu nalezy zmienic wartosc zmiennej TEST
// TEST=1 - siatka 4x4, TEST=2- siatka 4x4MIX, TEST=3, sitka 31x31
#define TEST 2

struct Jakobian{
    double J[2][2];
    double J_odwr[2][2];
    double detJ;

    void oblicz() {
        double a = J[0][0];
        double b = J[0][1];
        double c = J[1][0];
        double d = J[1][1];

        // Obliczanie wyznacznika
        detJ = a * d - b * c;

        // Sprawdzenie, czy wyznacznik jest różny od zera (macierz jest odwracalna)
        if (detJ != 0) {
            // Obliczanie odwrotności macierzy
            J_odwr[0][0] = d / detJ;
            J_odwr[0][1] = -b / detJ;
            J_odwr[1][0] = -c / detJ;
            J_odwr[1][1] = a / detJ;
        } else {
            std::cerr << "Macierz J nie jest odwracalna (detJ = 0)." << std::endl;
        }
    }

    // Funkcja do wypisywania danych macierzy
    void wypisz() const {
        std::cout << "\nMacierz J:" << std::endl;
        std::cout << J[0][0] << " " << J[0][1] << std::endl;
        std::cout << J[1][0] << " " << J[1][1] << std::endl;

        std::cout << "Wyznacznik detJ: " << detJ << std::endl;

        if (detJ != 0) {
            std::cout << "Macierz odwrotna J_odwr:" << std::endl;
            std::cout << J_odwr[0][0] << " " << J_odwr[0][1] << std::endl;
            std::cout << J_odwr[1][0] << " " << J_odwr[1][1] << std::endl;
        } else {
            std::cout << "Brak macierzy odwrotnej, ponieważ wyznacznik wynosi 0." << std::endl;
        }
    }
    void wypisz_odwr() const {
         if (detJ != 0){
            std::cout << "Macierz odwrotna J_odwr:" << std::endl;
            std::cout << J_odwr[0][0] << " " << J_odwr[0][1] << std::endl;
            std::cout << J_odwr[1][0] << " " << J_odwr[1][1] << std::endl;
        } else {
            std::cout << "Brak macierzy odwrotnej, ponieważ wyznacznik wynosi 0." << std::endl;
        }
    }
};
struct Macierz{
    double macierz[npc][4];

    void wypisz(std::string nazwa = "Macierz") const {
        std::cout << nazwa << ":" << std::endl;
        std::cout << std::fixed << std::setprecision(3);  // Ustawienie precyzji na 3 miejsca po przecinku
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                std::cout << std::setw(10) << macierz[i][j] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
};
#endif //HEADER_H
