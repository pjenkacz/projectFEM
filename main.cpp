#include <iostream>
#include <vector>
#include <fstream>
#include <math.h>
#include <sstream>
#include <iomanip>
#include "header.h"
//#define k 30
using namespace std;

#include <string>
#include <algorithm>

using namespace std;

struct Wezly{
    double X;
    double Y;
};

#if TEST == 1
std::string filename = "./testy/Test1_4_4.txt";
#endif
#if TEST == 2
std::string filename = "./testy/Test2_4_4_MixGrid.txt";
#endif
#if TEST == 3
std::string filename = "./testy/Test3_31_31_kwadrat.txt";
#endif

#if NPCBC == 2
double ksibc[NPCBC * 4] = {
    -1.0 / std::sqrt(3.0),  1.0 / std::sqrt(3.0),   // dolny bok (eta = -1)
     1.0,                   1.0,                   // prawa krawędź (ksi = 1)
     1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0),  // górny bok (eta = 1)
    -1.0,                  -1.0                    // lewy bok (ksi = -1)
};

double etabc[NPCBC * 4] = {
    -1.0, -1.0,                             // dolny bok (eta = -1, zmienna xi)
    -1.0 / std::sqrt(3.0),  1.0 / std::sqrt(3.0),   // prawa krawędź (ksi = 1, zmienna eta)
     1.0,  1.0,                               // górny bok (eta = 1, zmienna xi)
     1.0 / std::sqrt(3.0), -1.0 / std::sqrt(3.0)    // lewy bok (ksi = -1, zmienna eta)
};
double weightsbc[NPCBC] = { 1.0, 1.0 };
#elif NPCBC == 3
double ksibc[NPCBC * 4] = { -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0),
                              1.0, 1.0, 1.0,
                              std::sqrt(3.0 / 5.0), 0.0, -std::sqrt(3.0 / 5.0),
                              -1.0, -1.0, -1.0 };
double etabc[NPCBC * 4] = { -1.0, -1.0, -1.0,
                              -std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0),
                               1.0, 1.0, 1.0,
                               std::sqrt(3.0 / 5.0), 0.0, -std::sqrt(3.0 / 5.0) };
double weightsbc[NPCBC] = { 5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0 };
#elif NPCBC == 4
double ksibc[NPCBC * 4] = {
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),

     1.0, 1.0, 1.0, 1.0,

     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),

    -1.0, -1.0, -1.0, -1.0
};

double etabc[NPCBC * 4] = {
    -1.0, -1.0, -1.0, -1.0,

    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),

     1.0, 1.0, 1.0, 1.0,

     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))
};

double weightsbc[NPCBC] = {
    (18.0 - std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 - std::sqrt(30.0)) / 36.0
};
#else
#error "Unsupported value for npcbc. Please define npc as 2 or 3."
#endif

#if npc == 4
const int numPoints = 4;
Wezly wezly[numPoints] = {
    {-0.5773502692, -0.5773502692},
    { 0.5773502692, -0.5773502692},
    { 0.5773502692,  0.5773502692},
    {-0.5773502692,  0.5773502692}
};
Wezly wezlyB[numPoints] = {
    {-0.5773502692, -0.5773502692},
    { 0.5773502692, -0.5773502692},
    { 0.5773502692,  0.5773502692},
    {-0.5773502692,  0.5773502692}
};
double waga[3] = {1.0, 1.0,1.0};
#elif npc == 9
const int numPoints = 9;
Wezly wezly[numPoints] = {
    {-0.7745966692, -0.7745966692},
    {-0.7745966692, 0.0},
    {-0.7745966692, 0.7745966692},
    {0.0,  -0.7745966692     },
    {0.0,           0.0         },
    {0.0,  0.7745966692        },
    {0.7745966692,  -0.7745966692},
    {0.7745966692,           0.0},
    {0.7745966692,  0.7745966692}
};
double waga[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
#elif npc == 16
const int numPoints = 16;
double value[4] = {
    -std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
    -std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 - 2.0 / 7.0 * std::sqrt(6.0 / 5.0)),
     std::sqrt(3.0 / 7.0 + 2.0 / 7.0 * std::sqrt(6.0 / 5.0))
};

Wezly wezly[numPoints] = {
    { value[0], value[0] },
    { value[0], value[1] },
    { value[0], value[2] },
    { value[0], value[3] },

    { value[1], value[0] },
    { value[1], value[1] },
    { value[1], value[2] },
    { value[1], value[3] },

    { value[2], value[0] },
    { value[2], value[1] },
    { value[2], value[2] },
    { value[2], value[3] },

    { value[3], value[0] },
    { value[3], value[1] },
    { value[3], value[2] },
    { value[3], value[3] }
};

double waga[4] = {
    (18.0 - std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 + std::sqrt(30.0)) / 36.0,
    (18.0 - std::sqrt(30.0)) / 36.0
};


#else
#error "Unsupported value for npc. Please define npc as 4 or 9."
#endif
typedef vector<double> Vector;
typedef vector<Vector> Matrix;
void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name);
void printVector(const std::vector<double>& vec, const std::string& name);
Matrix createGlobal(const Matrix &Hglobal, const Matrix &Cglobal,  double dt);
Vector createRS(const Matrix &Cglobal, double dt,  const Vector &Tn, const Vector &Pglobal);
Vector gaussElimination(Matrix A, Vector b);
void solve(const Matrix &Hglobal,const Matrix &Cglobal,const Vector &Pglobal,double dt, double totalTime,double initialTemp);

struct ElemUniv{
    double dNdKsi[npc][4];
    double dNdNi[npc][4];

};
struct ElemUnivBC {
    double N[NPCBC * 4][4];  // 8 wierszy (każdy „punkt brzegowy”), 4 kolumny (N1..N4)
};
struct Sumy {
    double dXdKsi[4];
    double dYdNi[4];
    double dXdNi[4];
    double dYdKsi[4];
};
struct Macierz;
struct Jakobian;

struct Punkt {
    double X;
    double Y;
    bool BC; // czy węzeł jest brzegowy
};
struct GlobalData {
    double SimulationTime;
    double SimulationStepTime;
    double Conductivity;
    double Alfa;
    double Tot;
    double InitialTemp;
    double Density;
    double SpecificHeat;
    int nodesNumber;
    int elementsNumber;
};

struct Element {
    int ids[4];
    double H[4][4];
    double Hbc[4][4];
    double P[4];
    double C[4][4];
};
bool wczytajDane(const std::string& nazwaPliku, GlobalData& globalData,
                 std::vector<Punkt>& punkty, std::vector<Element>& elementy);


int main() {
    GlobalData globalData;
    vector<Punkt> punkty;
    vector<Element> elementy;
    //elementy.resize(globalData.elementsNumber);
    vector<int> BCnodes;
    // if(numPoints == 16){
    //     int idx = 0;
    //     for(int i=0; i<4; i++){
    //         for(int j=0; j<4; j++){
    //             wezly[idx].X = value[i];
    //             wezly[idx].Y = value[j];
    //             idx++;
    //         }
    //     }
    // }
    if (wczytajDane(filename, globalData, punkty, elementy)) {
        cout << "Dane wczytane poprawnie!" << endl;
        cout << "SimulationTime: " << globalData.SimulationTime << endl;
        cout << "Nodes: " << globalData.nodesNumber << ", Elements: " << globalData.elementsNumber << endl;
        cout << "Pierwszy punkt: X=" << punkty[0].X << ", Y=" << punkty[0].Y << ", BC=" << (punkty[0].BC ? "true" : "false") << endl;
        cout << endl;
    } else {
        cerr << "Błąd wczytywania danych." << endl;
    }

    cout << "Wczytano dane globalne:\n";
    cout << "k = " << globalData.Conductivity << "\n";
    //vector<Macierz> H_matrices;// Wektor, w którym przechowamy macierze wynikowe dla każdego elementu
    //vector<Macierz> Hbc_matrices;
    //H_matrices.resize(globalData.elementsNumber);
    //Hbc_matrices.resize(globalData.elementsNumber);
    std::vector<Macierz> Hbc_all(globalData.elementsNumber);
    std::vector<std::vector<double>> HGlobal;
    std::vector<std::vector<double>> CGlobal;
    std::vector<double> PGlobal;
    HGlobal.resize(globalData.nodesNumber, vector<double>(globalData.nodesNumber, 0.0));
    CGlobal.resize(globalData.nodesNumber, vector<double>(globalData.nodesNumber, 0.0));
    PGlobal.resize(globalData.nodesNumber, 0.0);
    cout<<"Wspolrzedzne wezlow : \n";
    for(auto x: wezly){
      cout<<x.X<<" "<<x.Y<<"\n";
    }

    //wypisanie wspolrzednych wczytanych punktow
    for (int i = 0; i < globalData.nodesNumber; ++i) {
        std::cout << "Punkt " << i + 1 << ": X = " << punkty[i].X
                  << ", Y = " << punkty[i].Y << '\n';
    }

    //wyzerowanie macierzy
    for (int e = 0; e < globalData.elementsNumber; e++) {
        // Wypełnienie struktury Macierz na 0
        for (int i=0; i<4; i++){
            for(int j=0;j<4;j++){
                Hbc_all[e].macierz[i][j]=0.0;
            }
        }

    }


    double dxdKSI, dydKSI, dxdNI, dydNI;
    ElemUniv pochodne; Sumy sumy; ElemUnivBC bcData;

    // pochodne sa takie same dla wszystkich elementow - nie zaleza od wspolrzednych punktow tylko wezlow
    //roznia sie zaleznie od 4/9/16 punktowego schematu
    //liczenie pochodnych "z palca"
    for (int i = 0; i < npc; i++) {
    double ksi = wezly[i].X;
    double ni = wezly[i].Y;

    pochodne.dNdKsi[i][0] = -0.25 * (1.0 - ni);
    pochodne.dNdKsi[i][1] =  0.25 * (1.0 - ni);
    pochodne.dNdKsi[i][2] =  0.25 * (1.0 + ni);
    pochodne.dNdKsi[i][3] = -0.25 * (1.0 + ni);

    pochodne.dNdNi[i][0] = -0.25 * (1.0 - ksi);
    pochodne.dNdNi[i][1] = -0.25 * (1.0 + ksi);
    pochodne.dNdNi[i][2] =  0.25 * (1.0 + ksi);
    pochodne.dNdNi[i][3] =  0.25 * (1.0 - ksi);
    }

    for (int i = 0; i < NPCBC * 4; ++i) {
        double ksi = ksibc[i];
        double eta = etabc[i];

        bcData.N[i][0] = 0.25 * (1.0 - ksi) * (1.0 - eta); // N1
        bcData.N[i][1] = 0.25 * (1.0 + ksi) * (1.0 - eta); // N2
        bcData.N[i][2] = 0.25 * (1.0 + ksi) * (1.0 + eta); // N3
        bcData.N[i][3] = 0.25 * (1.0 - ksi) * (1.0 + eta); // N4
    }

    //od tego momentu konieczne jest osobne liczenie dla kazdego elementu

    Jakobian Jakobiany[npc];
    for (int e = 0; e < globalData.elementsNumber; e++) {
    // Pobieramy id węzłów elementu
    int n1 = elementy[e].ids[0] - 1;
    int n2 = elementy[e].ids[1] - 1;
    int n3 = elementy[e].ids[2] - 1;
    int n4 = elementy[e].ids[3] - 1;

    Punkt localPunkty[4] = { punkty[n1], punkty[n2], punkty[n3], punkty[n4] };

    // Obliczamy Jakobiany dla tego elementu
    //liczenie sum po dx dy "z palca"
    for (int j = 0; j < npc; j++) {
        double dxdKSI = pochodne.dNdKsi[j][0]*localPunkty[0].X + pochodne.dNdKsi[j][1]*localPunkty[1].X +
                        pochodne.dNdKsi[j][2]*localPunkty[2].X + pochodne.dNdKsi[j][3]*localPunkty[3].X;

        double dydKSI = pochodne.dNdKsi[j][0]*localPunkty[0].Y + pochodne.dNdKsi[j][1]*localPunkty[1].Y +
                        pochodne.dNdKsi[j][2]*localPunkty[2].Y + pochodne.dNdKsi[j][3]*localPunkty[3].Y;

        double dxdNI = pochodne.dNdNi[j][0]*localPunkty[0].X + pochodne.dNdNi[j][1]*localPunkty[1].X +
                       pochodne.dNdNi[j][2]*localPunkty[2].X + pochodne.dNdNi[j][3]*localPunkty[3].X;

        double dydNI = pochodne.dNdNi[j][0]*localPunkty[0].Y + pochodne.dNdNi[j][1]*localPunkty[1].Y +
                       pochodne.dNdNi[j][2]*localPunkty[2].Y + pochodne.dNdNi[j][3]*localPunkty[3].Y;

        Jakobiany[j].J[0][0] = dxdKSI;
        Jakobiany[j].J[1][0] = dxdNI;
        Jakobiany[j].J[0][1] = dydKSI;
        Jakobiany[j].J[1][1] = dydNI;
        Jakobiany[j].oblicz();
    }

    // Obliczamy pochodne dN/dx i dN/dy dla punktów całkowania tego elementu
    double dNdx[npc][4];
    double dNdy[npc][4];

    for (int i = 0; i < npc; ++i) {
        for (int j = 0; j < 4; ++j) {
            dNdx[i][j] = Jakobiany[i].J_odwr[0][0]*pochodne.dNdKsi[i][j] + Jakobiany[i].J_odwr[0][1]*pochodne.dNdNi[i][j];
            dNdy[i][j] = Jakobiany[i].J_odwr[1][0]*pochodne.dNdKsi[i][j] + Jakobiany[i].J_odwr[1][1]*pochodne.dNdNi[i][j];
        }
    }

    // Liczymy macierze H dla tego elementu
    Macierz H_X[npc] = {0};
    Macierz H_Y[npc] = {0};
    Macierz H_Suma[npc] = {0};
    Macierz C_Suma[npc] = {0};
    for (int s = 0; s < npc; s++) {

        double N[4];
    {
        double ksi = wezly[s].X;
        double eta = wezly[s].Y;
        N[0] = 0.25 * (1.0 - ksi) * (1.0 - eta);
        N[1] = 0.25 * (1.0 + ksi) * (1.0 - eta);
        N[2] = 0.25 * (1.0 + ksi) * (1.0 + eta);
        N[3] = 0.25 * (1.0 - ksi) * (1.0 + eta);
    }
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                H_X[s].macierz[i][j] = dNdx[s][i] * dNdx[s][j] * Jakobiany[s].detJ;
                H_Y[s].macierz[i][j] = dNdy[s][i] * dNdy[s][j] * Jakobiany[s].detJ;
                H_Suma[s].macierz[i][j] = (H_X[s].macierz[i][j] + H_Y[s].macierz[i][j]) * globalData.Conductivity;
                C_Suma[s].macierz[i][j] = N[i] * N[j] * Jakobiany[s].detJ * globalData.Density * globalData.SpecificHeat;
            }
        }
    }

    Macierz wynik = {0};
    Macierz wynikC = {0};
    double w1, w2;
    if (npc == 4) {
        for (int s = 0; s < npc; ++s) {
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    wynik.macierz[i][j] += H_Suma[s].macierz[i][j]; // wagi = 1 * 1
                    wynikC.macierz[i][j] += C_Suma[s].macierz[i][j];
                }
            }
        }
    } else if (npc == 9) {
        for (int s = 0; s < npc; ++s) {
            int row = s / 3; // wiersz w macierzy 3x3
            int col = s % 3; // kolumna w macierzy 3x3

            w1 = waga[row];
            w2 = waga[col];

            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    wynik.macierz[i][j] += H_Suma[s].macierz[i][j] * w1 * w2;
                    wynikC.macierz[i][j] += C_Suma[s].macierz[i][j] * w1 * w2;
                }
            }
        }
    }
    else if (npc == 16) {
        for (int s = 0; s < npc; ++s) {
            int row = s / 4;
            int col = s % 4;
            w1 = waga[row];
            w2 = waga[col];
            for (int i = 0; i < 4; ++i) {
                for (int j = 0; j < 4; ++j) {
                    wynik.macierz[i][j]  += H_Suma[s].macierz[i][j] * w1 * w2;
                    wynikC.macierz[i][j] += C_Suma[s].macierz[i][j] * w1 * w2;
                }
            }
        }
    }
    // Zapisujemy obliczone macierze H i C do struktury elementu
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            elementy[e].H[i][j] = wynik.macierz[i][j];
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
           elementy[e].C[i][j] = wynikC.macierz[i][j];
    }
    }

    // Obliczanie Macierzy Hbc oraz wektora P
    Macierz Hbc_local = {0};
    double P_local[4] = {0.0, 0.0, 0.0, 0.0};
    double N[4];
    for (int edge = 0; edge < 4; ++edge) {

        Punkt p1 = localPunkty[edge];
        Punkt p2 = localPunkty[(edge + 1) % 4];

        // Długość boku liczona z Pitagorasa
        double dx = p2.X - p1.X;
        double dy = p2.Y - p1.Y;
        double length = std::sqrt(dx*dx + dy*dy);
        double detJ_1D = length / 2.0;
        double w_i   = 1.0;
        if(NPCBC == 2){
        if(p1.BC && p2.BC){

        for (int iBC = 0; iBC < NPCBC; ++iBC) {
            int idxBC = edge * NPCBC + iBC;

            for(int iN=0; iN<4; iN++){
                N[iN] = bcData.N[idxBC][iN];
            }
            for (int a = 0; a < 4; ++a) {
                for (int b = 0; b < 4; ++b) {
                    Hbc_local.macierz[a][b] += globalData.Alfa * N[a] * N[b] * detJ_1D * w_i;
                        }
                }
             for (int a = 0; a < 4; ++a) {
                 P_local[a] += globalData.Alfa * globalData.Tot * N[a] * detJ_1D * w_i;
                }
        } } }
        else if(NPCBC == 3){
            if(p1.BC && p2.BC){

            for (int iBC = 0; iBC < NPCBC; ++iBC) {
                    int idxBC = edge * NPCBC + iBC;
                    double ksi  = ksibc[idxBC];
                    double eta  = etabc[idxBC];
                    double w    = weightsbc[iBC];    // wagi

                   // double N[4];
                    N[0] = 0.25 * (1 - ksi) * (1 - eta);
                    N[1] = 0.25 * (1 + ksi) * (1 - eta);
                    N[2] = 0.25 * (1 + ksi) * (1 + eta);
                    N[3] = 0.25 * (1 - ksi) * (1 + eta);
            for (int a = 0; a < 4; ++a) {
                for (int b = 0; b < 4; ++b) {
                    Hbc_local.macierz[a][b] += globalData.Alfa * N[a] * N[b] * detJ_1D * w;
                        }
                    }
            for (int a = 0; a < 4; a++) {
            P_local[a] += globalData.Alfa * globalData.Tot * N[a] * detJ_1D * w;
                }
        } } }
        else if(NPCBC == 4){
        if(p1.BC && p2.BC){
        for (int iBC = 0; iBC < NPCBC; ++iBC) {
            int idxBC = edge * NPCBC + iBC;
            double ksi  = ksibc[idxBC];
            double eta  = etabc[idxBC];
            double w    = weightsbc[iBC];
            double N[4];
                    N[0] = 0.25 * (1 - ksi) * (1 - eta);
                    N[1] = 0.25 * (1 + ksi) * (1 - eta);
                    N[2] = 0.25 * (1 + ksi) * (1 + eta);
                    N[3] = 0.25 * (1 - ksi) * (1 + eta);
            for (int a = 0; a < 4; ++a) {
                for (int b = 0; b < 4; ++b) {
                    Hbc_local.macierz[a][b] += globalData.Alfa * N[a] * N[b] * detJ_1D * w;
                        }
                    }
            for (int a = 0; a < 4; a++) {
            P_local[a] += globalData.Alfa * globalData.Tot * N[a] * detJ_1D * w;
                }
        }
    }
}

    }
    // Zapisujemy macierz Hbc oraz wektor P do struktury lementu
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            elementy[e].Hbc[i][j] = Hbc_local.macierz[i][j];
        }
        elementy[e].P[i] = P_local[i];
    }
    //Wypisanie macierzy H dla każdego elementu
    std::cout << "Macierz H dla elementu " << e+1 << ":" << std::endl;
    wynik.wypisz("Macierz H elementu");
    std::cout << "Macierz Hbc dla elementu " << e+1 << ":" << std::endl;
    Hbc_local.wypisz("Macierz Hbc elementu");
    std::cout << "Macierz C dla elementu " << e+1 << ":" << std::endl;
    wynikC.wypisz("Macierz C elementu");
    }
    //Agregacja macierzy dla kazdego elementu do macierzy globalnych

    for (int i = 0; i < globalData.elementsNumber; ++i) {

    for (int j = 0; j < 4; ++j) {
        int globalJ = elementy[i].ids[j] - 1;  // -1 bo numeracja w pliku zaczyna sie od 1
        for (int k = 0; k < 4; ++k) {
            int globalK = elementy[i].ids[k] - 1;

            // Dodajemy macierze elementow H do macierzy globalnej H
            HGlobal[globalJ][globalK] += elementy[i].H[j][k];
            // Dodajemy macierze elementow Hbc do macierzy globalnej H
            HGlobal[globalJ][globalK] += elementy[i].Hbc[j][k];
            // Dodajemy macierze elementow C do macierzy globalnej H
            CGlobal[globalJ][globalK] += elementy[i].C[j][k];
        }
        // Dodajemy do wektora globalnego
        PGlobal[globalJ] += elementy[i].P[j];
    }
}
    // printMatrix(HGlobal,"HGlobal");
    // printMatrix(CGlobal,"CGlobal");
    // printVector(PGlobal, "PGlobal");

    solve(HGlobal, CGlobal, PGlobal, globalData.SimulationStepTime, globalData.SimulationTime, globalData.InitialTemp);
    return 0;
}



bool wczytajDane(const std::string& nazwaPliku, GlobalData& globalData,
                 std::vector<Punkt>& punkty, std::vector<Element>& elementy) {
    std::ifstream plik(nazwaPliku);
    if (!plik.is_open()) {
        std::cerr << "Nie można otworzyć pliku: " << nazwaPliku << std::endl;
        return false;
    }

    std::string linia;
    bool nodesSection = false;
    bool elementsSection = false;

    while (std::getline(plik, linia)) {
        if (linia.empty()) continue;

        if (linia.find("*Node") != std::string::npos) {
            nodesSection = true;
            break;
        }

        std::istringstream iss(linia);
        std::string paramName;
        iss >> paramName;
        if (paramName == "SimulationTime") {
            iss >> globalData.SimulationTime;
        } else if (paramName == "SimulationStepTime") {
            iss >> globalData.SimulationStepTime;
        } else if (paramName == "Conductivity") {
            iss >> globalData.Conductivity;
        } else if (paramName == "Alfa") {
            iss >> globalData.Alfa;
        } else if (paramName == "Tot") {
            iss >> globalData.Tot;
        } else if (paramName == "InitialTemp") {
            iss >> globalData.InitialTemp;
        } else if (paramName == "Density") {
            iss >> globalData.Density;
        } else if (paramName == "SpecificHeat") {
            iss >> globalData.SpecificHeat;
        } else if (paramName == "Nodes") {
            std::string tmp;
            iss >> tmp >> globalData.nodesNumber;
        } else if (paramName == "Elements") {
            std::string tmp;
            iss >> tmp >> globalData.elementsNumber;
        }
    }

    if (!nodesSection) {
        std::cerr << "Nie znaleziono sekcji *Node w pliku!" << std::endl;
        return false;
    }

    // Wczytanie węzłów
    punkty.resize(globalData.nodesNumber);
    for (int i = 0; i < globalData.nodesNumber; i++) {
        if (!std::getline(plik, linia)) {
            std::cerr << "Oczekiwano więcej węzłów." << std::endl;
            return false;
        }
        std::replace(linia.begin(), linia.end(), ',', ' ');
        std::istringstream iss(linia);
        int id;
        Punkt p;
        iss >> id >> p.X >> p.Y;
        p.BC = false;
        punkty[id - 1] = p;
    }

    // Wczytanie elementów
    if (!std::getline(plik, linia) || linia.find("*Element") == std::string::npos) {
        std::cerr << "Nie odnaleziono linii *Element, type=DC2D4" << std::endl;
        return false;
    }

    elementy.resize(globalData.elementsNumber);
    for (int i = 0; i < globalData.elementsNumber; i++) {
        if (!std::getline(plik, linia)) {
            std::cerr << "Oczekiwano więcej elementów." << std::endl;
            return false;
        }
        std::replace(linia.begin(), linia.end(), ',', ' ');
        std::istringstream iss(linia);
        int elemID;
        Element elem;
        iss >> elemID >> elem.ids[0] >> elem.ids[1] >> elem.ids[2] >> elem.ids[3];
        elementy[i] = elem;
    }

    // Odczyt *BC
    if (!std::getline(plik, linia) || linia.find("*BC") == std::string::npos) {
        std::cerr << "Nie znaleziono sekcji *BC" << std::endl;
        return false;
    }

    // Kolejna linia zawiera węzły brzegowe
    if (!std::getline(plik, linia)) {
        std::cerr << "Nie podano węzłów brzegowych." << std::endl;
        return false;
    }
    std::replace(linia.begin(), linia.end(), ',', ' ');
    {
        std::istringstream iss(linia);
        int bcNode;
        while (iss >> bcNode) {
            if (bcNode >=1 && bcNode <= globalData.nodesNumber) {
                punkty[bcNode - 1].BC = true;
            }
        }
    }

    return true;
}
void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name) {
    std::cout << name << ":" << std::endl;
    for (const auto& row : matrix) {
        for (const auto& value : row) {
            std::cout << std::setw(4) << std::fixed << std::setprecision(3) << value << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
void printVector(const std::vector<double>& vec, const std::string& name) {
    std::cout << name << ":" << std::endl;
    for (const auto& value : vec) {
        std::cout << std::setw(4) << std::fixed << std::setprecision(3) << value << " ";
    }
    std::cout << std::endl << std::endl;
}

Matrix createGlobal(const Matrix &Hglobal,
                    const Matrix &Cglobal,
                    double dt)
{
    int N = (int)Hglobal.size();
    Matrix Global(N, Vector(N, 0.0));
    // Tworzymy macierz Globalna wg wzoru
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Global[i][j] = Hglobal[i][j] + Cglobal[i][j] / dt;
        }
    }
    return Global;
}

Vector createRS(const Matrix &Cglobal, double dt,const Vector &Tn,const Vector &Pglobal)
{
    int N = (int)Cglobal.size();
    Vector RS(N, 0.0);

    // Tworzymy wektor RS wg wzoru
    for (int i = 0; i < N; i++) {
        double sumC = 0.0;
        for (int j = 0; j < N; j++) {
            sumC += (Cglobal[i][j] / dt) * Tn[j];
        }
        RS[i] = sumC + Pglobal[i];
    }
    return RS;
}
Vector gaussElimination(Matrix A, Vector b) {
    int N = (int)A.size();

    //  Eliminacja przód macierz trojkatna gorna
    for (int i = 0; i < N; i++) {
        // Dzielenie wiersza i (żeby A[i][i] było =1)
        double diag = A[i][i];
        for (int col = i; col < N; col++) {
            A[i][col] /= diag;
        }
        b[i] /= diag;

        // Odejmowanie wierszy
        for (int row = i + 1; row < N; row++) {
            double factor = A[row][i];
            for (int col = i; col < N; col++) {
                A[row][col] -= factor * A[i][col];
            }
            b[row] -= factor * b[i];
        }
    }

    //  Eliminacja wstecz
    Vector x(N, 0.0);
    for (int i = N - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < N; j++) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum;
    }
    return x;
}
void solve(const Matrix &Hglobal,
                      const Matrix &Cglobal,
                      const Vector &Pglobal,
                      double dt,
                      double totalTime,
                      double initialTemp)
{
    int N = (int)Hglobal.size();
    Matrix Global = createGlobal(Hglobal, Cglobal, dt); // Obliczamy macierz Global = H + C/dt
    Vector T0(N, initialTemp);                          // Wektor temperatur
    Vector T(N, 0.0);                                    //wektor pomocniczy
    Vector RS = createRS(Cglobal, dt, T0, Pglobal);  // Obliczamy wektor RS = (C/dt)*T0 + Pglobal - prawa strona rownania
    //  printVector(RS, "Iteracja 0");
    //  printMatrix(Global, "Iteracja 0");
    for (double currentTime = dt; currentTime <= totalTime; currentTime += dt) {
        Vector RS = createRS(Cglobal, dt, T0, Pglobal);  // Obliczamy wektor RS = (C/dt)*T0 + Pglobal - prawa strona rownania
        T = gaussElimination(Global, RS); // Rozwiązujemy układ: Global * T = RS

        double minT = *std::min_element(T.begin(), T.end());
        double maxT = *std::max_element(T.begin(), T.end());
        std::cout << std::fixed << std::setprecision(4);
        cout << "Czas = " << currentTime
             << " [s], Tmin = " << minT
             << ", Tmax = " << maxT << endl;

        // przygotowanie do następnej iteracji
        T0 = T;
    }
}
