#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

const int N = 5;
const double EPSILON = 1e-12;

void dekompozycja_LU(double A[N][N], int wektor_permutacji[N]) {
    for (int i = 0; i < N; i++)
        wektor_permutacji[i] = i;

    for (int kolumna = 0; kolumna < N; kolumna++) {
        double max_val = fabs(A[wektor_permutacji[kolumna]][kolumna]);
        int wiersz_max = kolumna;

        // Szukamy wiersza, który ma największą wartość bezwzględną w bieżącej kolumnie
        for (int wiersz = kolumna + 1; wiersz < N; wiersz++) {
            double aktualna_wartosc = fabs(A[wektor_permutacji[wiersz]][kolumna]);
            if (aktualna_wartosc > max_val) {
                max_val = aktualna_wartosc;
                wiersz_max = wiersz;
            }
        }

        if (max_val < EPSILON) {
            cout << "Macierz osobliwa lub bliska osobliwej!" << endl;
            return;
        }

        // Zamiana wierszy w wektorze permutacji
        swap(wektor_permutacji[kolumna], wektor_permutacji[wiersz_max]);

        double element_glowny = A[wektor_permutacji[kolumna]][kolumna];

        for (int wiersz = kolumna + 1; wiersz < N; wiersz++) {
            double mnoznik = A[wektor_permutacji[wiersz]][kolumna] / element_glowny;
            A[wektor_permutacji[wiersz]][kolumna] = mnoznik;

            for (int j = kolumna + 1; j < N; j++)
                A[wektor_permutacji[wiersz]][j] -= mnoznik * A[wektor_permutacji[kolumna]][j];
        }
    }
}

void rozwiaz_uklad(double A[N][N], int wektor_permutacji[N], double b[N], double x[N]) {
    double y[N];

    // Podstawianie w przód (Ly = b)
    for (int i = 0; i < N; i++) {
        y[i] = b[wektor_permutacji[i]];
        for (int j = 0; j < i; j++)
            y[i] -= A[wektor_permutacji[i]][j] * y[j];
    }

    // Podstawianie wstecz (Ux = y)
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++)
            x[i] -= A[wektor_permutacji[i]][j] * x[j];
        x[i] /= A[wektor_permutacji[i]][i];
    }
}

int main() {
    // Macierz A
    double A[N][N] = {
        { 5, 4, 3, 2, 1},
        {10, 8, 7, 6, 5},
        {-1, 2, -3, 4, -5},
        { 6, 5, -4, 3, -2},
        { 1, 2, 3, 4, 5}
    };

    // Wektor b
    double b[N] = {37, 99, -9, 12, 53};
    int wektor_permutacji[N];
    double x[N];

    dekompozycja_LU(A, wektor_permutacji);
    rozwiaz_uklad(A, wektor_permutacji, b, x);

    cout << "Rozwiazanie ukladu (x):" << endl;
    for (int i = 0; i < N; i++)
        cout << "x[" << i + 1 << "] = " << x[i] << endl;

    return 0;
}