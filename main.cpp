#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <string>

using namespace std;

const int ROZMIAR = 5; // liczba rownan / niewiadomych

// Macierz wspolczynnikow A i wektor prawych stron b
double macierzA[ROZMIAR][ROZMIAR] = {
    {50,  5,  4,  3,  2},
    { 1, 40,  1,  2,  3},
    { 4,  5, 30, -5, -4},
    {-3, -2, -1, 20,  0},
    { 1,  2,  3,  4, 30}
};

double prawaStrona[ROZMIAR] = {140, 67, 62, 89, 153};

double przyblizenieStartowe[ROZMIAR] = {6, 6, 6, 6, 6};

// Oblicza norme L2 residuum: ||b - A*x||
double normaResiduum(double rozwiazanie[ROZMIAR]) {
    double suma = 0.0;
    for (int wiersz = 0; wiersz < ROZMIAR; wiersz++) {
        double skladnik_residuum = prawaStrona[wiersz];
        for (int kolumna = 0; kolumna < ROZMIAR; kolumna++)
            skladnik_residuum -= macierzA[wiersz][kolumna] * rozwiazanie[kolumna];
        suma += skladnik_residuum * skladnik_residuum;
    }
    return sqrt(suma);
}

// Oblicza norme L2 roznicy dwoch wektorow: ||x_nowe - x_poprzednie|| (estymator bledu)
double normaRoznicy(double rozwiazanie[ROZMIAR], double poprzednie[ROZMIAR]) {
    double suma = 0.0;
    for (int i = 0; i < ROZMIAR; i++) {
        double roznica = rozwiazanie[i] - poprzednie[i];
        suma += roznica * roznica;
    }
    return sqrt(suma);
}

// Wypisuje wyniki biezacej iteracji
void wypiszIteracje(int nrIteracji, double rozwiazanie[ROZMIAR], double estymatorBledu, double residuum) {
    cout << "  Iter " << setw(4) << nrIteracji << " | x = [";
    for (int i = 0; i < ROZMIAR; i++)
        cout << setw(10) << fixed << setprecision(6) << rozwiazanie[i] << (i < ROZMIAR-1 ? ", " : "");
    cout << "] | ||dx|| = " << scientific << setprecision(4) << estymatorBledu
         << " | ||r|| = "  << scientific << setprecision(4) << residuum << "\n";
}

// -------------------------
// METODA JACOBIEGO
// -------------------------
void jacobi(double tolerancja1, double tolerancja2, int maksIteracji) {
    cout << "\n========== METODA JACOBIEGO ==========\n";
    cout << "Kryteria stopu: (1) ||dx||<" << tolerancja1
         << "  (2) ||r||<" << tolerancja2
         << "  (3) iter>" << maksIteracji << "\n\n";

    double rozwiazanie[ROZMIAR], noweRozwiazanie[ROZMIAR];
    for (int i = 0; i < ROZMIAR; i++) rozwiazanie[i] = przyblizenieStartowe[i];

    int nrIteracji        = 0;
    double estymatorBledu = 1e12;
    double residuum       = 1e12;

    while (nrIteracji < maksIteracji) {
        // Wyznaczamy nowe przyblizenie (uzywamy TYLKO starych wartosci - cecha Jacobiego)
        for (int wiersz = 0; wiersz < ROZMIAR; wiersz++) {
            double sumaOffDiag = prawaStrona[wiersz];
            for (int kolumna = 0; kolumna < ROZMIAR; kolumna++)
                if (kolumna != wiersz)
                    sumaOffDiag -= macierzA[wiersz][kolumna] * rozwiazanie[kolumna];
            noweRozwiazanie[wiersz] = sumaOffDiag / macierzA[wiersz][wiersz];
        }

        estymatorBledu = normaRoznicy(noweRozwiazanie, rozwiazanie);
        residuum       = normaResiduum(noweRozwiazanie);
        nrIteracji++;

        for (int i = 0; i < ROZMIAR; i++) rozwiazanie[i] = noweRozwiazanie[i];

        wypiszIteracje(nrIteracji, rozwiazanie, estymatorBledu, residuum);

        // Sprawdzenie kryteriow stopu
        if (estymatorBledu < tolerancja1) {
            cout << "  >> Zbieznosc: ||dx|| < " << tolerancja1 << " po " << nrIteracji << " iteracjach.\n";
            break;
        }
        if (residuum < tolerancja2) {
            cout << "  >> Zbieznosc: ||r|| < " << tolerancja2 << " po " << nrIteracji << " iteracjach.\n";
            break;
        }
        if (nrIteracji >= maksIteracji) {
            cout << "  >> Osiagnieto maks. liczbe iteracji: " << maksIteracji << "\n";
        }
    }
    cout << "  Koncowe rozwiazanie: [";
    for (int i = 0; i < ROZMIAR; i++)
        cout << fixed << setprecision(6) << rozwiazanie[i] << (i < ROZMIAR-1 ? ", " : "");
    cout << "]\n";
    cout << "  Liczba iteracji: " << nrIteracji << "\n";
}

// -------------------------
// METODA GAUSSA-SEIDELA
// -------------------------
void gaussSeidel(double tolerancja1, double tolerancja2, int maksIteracji) {
    cout << "\n========== METODA GAUSSA-SEIDELA ==========\n";
    cout << "Kryteria stopu: (1) ||dx||<" << tolerancja1
         << "  (2) ||r||<" << tolerancja2
         << "  (3) iter>" << maksIteracji << "\n\n";

    double rozwiazanie[ROZMIAR], poprzednieRozwiazanie[ROZMIAR];
    for (int i = 0; i < ROZMIAR; i++) rozwiazanie[i] = przyblizenieStartowe[i];

    int nrIteracji        = 0;
    double estymatorBledu = 1e12;
    double residuum       = 1e12;

    while (nrIteracji < maksIteracji) {
        for (int i = 0; i < ROZMIAR; i++) poprzednieRozwiazanie[i] = rozwiazanie[i];

        // Aktualizujemy kazda skladowa od razu uzywajac swiezych wartosci (cecha Gaussa-Seidela)
        for (int wiersz = 0; wiersz < ROZMIAR; wiersz++) {
            double sumaOffDiag = prawaStrona[wiersz];
            for (int kolumna = 0; kolumna < ROZMIAR; kolumna++)
                if (kolumna != wiersz)
                    sumaOffDiag -= macierzA[wiersz][kolumna] * rozwiazanie[kolumna];
            rozwiazanie[wiersz] = sumaOffDiag / macierzA[wiersz][wiersz];
        }

        estymatorBledu = normaRoznicy(rozwiazanie, poprzednieRozwiazanie);
        residuum       = normaResiduum(rozwiazanie);
        nrIteracji++;

        wypiszIteracje(nrIteracji, rozwiazanie, estymatorBledu, residuum);

        if (estymatorBledu < tolerancja1) {
            cout << "  >> Zbieznosc: ||dx|| < " << tolerancja1 << " po " << nrIteracji << " iteracjach.\n";
            break;
        }
        if (residuum < tolerancja2) {
            cout << "  >> Zbieznosc: ||r|| < " << tolerancja2 << " po " << nrIteracji << " iteracjach.\n";
            break;
        }
        if (nrIteracji >= maksIteracji) {
            cout << "  >> Osiagnieto maks. liczbe iteracji: " << maksIteracji << "\n";
        }
    }
    cout << "  Koncowe rozwiazanie: [";
    for (int i = 0; i < ROZMIAR; i++)
        cout << fixed << setprecision(6) << rozwiazanie[i] << (i < ROZMIAR-1 ? ", " : "");
    cout << "]\n";
    cout << "  Liczba iteracji: " << nrIteracji << "\n";
}

// -------------------------
// METODA SOR (sukcesywna nadrelaksacja)
// -------------------------
void sor(double parametrOmega, double tolerancja1, double tolerancja2, int maksIteracji) {
    cout << "\n========== METODA SOR (omega = " << parametrOmega << ") ==========\n";
    cout << "Kryteria stopu: (1) ||dx||<" << tolerancja1
         << "  (2) ||r||<" << tolerancja2
         << "  (3) iter>" << maksIteracji << "\n\n";

    double rozwiazanie[ROZMIAR], poprzednieRozwiazanie[ROZMIAR];
    for (int i = 0; i < ROZMIAR; i++) rozwiazanie[i] = przyblizenieStartowe[i];

    int nrIteracji        = 0;
    double estymatorBledu = 1e12;
    double residuum       = 1e12;

    while (nrIteracji < maksIteracji) {
        for (int i = 0; i < ROZMIAR; i++) poprzednieRozwiazanie[i] = rozwiazanie[i];

        for (int wiersz = 0; wiersz < ROZMIAR; wiersz++) {
            double sumaOffDiag = prawaStrona[wiersz];
            for (int kolumna = 0; kolumna < ROZMIAR; kolumna++)
                if (kolumna != wiersz)
                    sumaOffDiag -= macierzA[wiersz][kolumna] * rozwiazanie[kolumna];
            // Krok Gaussa-Seidela
            double krokGS = sumaOffDiag / macierzA[wiersz][wiersz];
            // Relaksacja: srednia wazona miedzy starym przyblizeniem a krokiem GS
            rozwiazanie[wiersz] = (1.0 - parametrOmega) * poprzednieRozwiazanie[wiersz]
                                 + parametrOmega * krokGS;
        }

        estymatorBledu = normaRoznicy(rozwiazanie, poprzednieRozwiazanie);
        residuum       = normaResiduum(rozwiazanie);
        nrIteracji++;

        wypiszIteracje(nrIteracji, rozwiazanie, estymatorBledu, residuum);

        if (estymatorBledu < tolerancja1) {
            cout << "  >> Zbieznosc: ||dx|| < " << tolerancja1 << " po " << nrIteracji << " iteracjach.\n";
            break;
        }
        if (residuum < tolerancja2) {
            cout << "  >> Zbieznosc: ||r|| < " << tolerancja2 << " po " << nrIteracji << " iteracjach.\n";
            break;
        }
        if (nrIteracji >= maksIteracji) {
            cout << "  >> Osiagnieto maks. liczbe iteracji: " << maksIteracji << "\n";
        }
    }
    cout << "  Koncowe rozwiazanie: [";
    for (int i = 0; i < ROZMIAR; i++)
        cout << fixed << setprecision(6) << rozwiazanie[i] << (i < ROZMIAR-1 ? ", " : "");
    cout << "]\n";
    cout << "  Liczba iteracji: " << nrIteracji << "\n";
}

int main() {
    cout << "====================================================\n";
    cout << "  Uklad rownan liniowych Ax = b\n";
    cout << "  Rozmiar ukladu: " << ROZMIAR << "x" << ROZMIAR << "\n";
    cout << "  Przyblizenie startowe x0 = [6,6,6,6,6]\n";
    cout << "  Trzy kryteria stopu:\n";
    cout << "    (1) estymator bledu ||x^(n) - x^(n-1)|| < tol = 1e-6\n";
    cout << "    (2) norma residuum  ||b - Ax||           < tol = 1e-6\n";
    cout << "    (3) liczba iteracji > maksIteracji = 1000\n";
    cout << "====================================================\n";

    double tolerancja   = 1e-12;
    int    maksIteracji = 1000;

    jacobi(tolerancja, tolerancja, maksIteracji);
    gaussSeidel(tolerancja, tolerancja, maksIteracji);
    sor(0.5, tolerancja, tolerancja, maksIteracji);

    return 0;
}