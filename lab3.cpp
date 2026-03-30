#include <iostream>
#include <cmath>

using namespace std;

// Rownanie A: tanh(x) + 2x - 2 = 0
double rownanie_a(double t) {
    return tanh(t) + 2.0 * t - 2.0;
}

double pochodna_a(double t) {
    double ch = cosh(t);
    return 1.0 / (ch * ch) + 2.0;
}

double iteracja_a(double t) {
    return 1.0 - 0.5 * tanh(t);
}

double iteracja_prim_a(double t) {
    double ch = cosh(t);
    return -0.5 / (ch * ch);
}

// Rownanie B: sinh(x) + x/4 - 1 = 0
double rownanie_b(double t) {
    return sinh(t) + t / 4.0 - 1.0;
}

double pochodna_b(double t) {
    return cosh(t) + 0.25;
}

double iteracja_b(double t) {
    return asinh(1.0 - t / 4.0);
}

double iteracja_prim_b(double t) {
    double v = 1.0 - t / 4.0;
    return -0.25 / sqrt(v * v + 1.0);
}


double oblicz_f(char p, double t) {
    return (p == 'a') ? rownanie_a(t) : rownanie_b(t);
}

double oblicz_fi(char p, double t) {
    return (p == 'a') ? iteracja_a(t) : iteracja_b(t);
}

double oblicz_fi_prim(char p, double t) {
    return (p == 'a') ? iteracja_prim_a(t) : iteracja_prim_b(t);
}

double oblicz_f_prim(char p, double t) {
    return (p == 'a') ? pochodna_a(t) : pochodna_b(t);
}


void metoda_picarda(char p, int limit, double tolx, double tolf) {
    double biezacy = (p == 'a') ? 0.0 : 0.5;

    for (int iter = 0; iter < limit; iter++) {
        double nastepny  = oblicz_fi(p, biezacy);
        double blad      = fabs(nastepny - biezacy);
        double res       = fabs(oblicz_f(p, nastepny));
        double wspQ      = fabs(oblicz_fi_prim(p, nastepny));

        printf("I: %d xn: %.16f xn+1: %.16f en: %.16f f(xn) %.16f %s\n",
               iter + 1, biezacy, nastepny, blad, res,
               wspQ < 1.0 ? "zbiezny" : "rozbiezny");

        if (blad <= tolx && res <= tolf) {
            printf("Oczekiwana wartosc uzyskana, konczymy iteracje. Wynik to: %.16f \n", nastepny);
            return;
        }

        biezacy = nastepny;
    }

    printf("Pierwiastek nie osiagnieto po %d iteracji, ostatnie przyblizenie to: %.16f \n", limit, biezacy);
}


void metoda_bisekcji(char p, int limit, double tolx, double tolf) {
    double lewy = 0.0, prawy = 1.0;

    double fl = oblicz_f(p, lewy);
    double fp = oblicz_f(p, prawy);

    if ((fl > 0.0 && fp > 0.0) || (fl < 0.0 && fp < 0.0)) {
        cout << "Funckja ta nie ma roznego znaku na koncach przedzialu \n" << endl;
        return;
    }

    double srodek;

    for (int iter = 0; iter < limit; iter++) {
        srodek          = (lewy + prawy) / 2.0;
        double blad     = (prawy - lewy) / 2.0;
        double res      = fabs(oblicz_f(p, srodek));

        printf("I: %d [a,b]: [%.14f,%.14f] xn: %.14f en: %.14f f(xn) %.14f \n",
               iter + 1, lewy, prawy, srodek, blad, res);

        if (blad <= tolx && res <= tolf) {
            printf("Oczekiwana wartosc uzyskana, konczymy iteracje. Wynik to: %.16f \n", srodek);
            return;
        }

        double fl_new = oblicz_f(p, lewy);
        double fs     = oblicz_f(p, srodek);

        if ((fl_new < 0.0 && fs > 0.0) || (fl_new > 0.0 && fs < 0.0))
            prawy = srodek;
        else
            lewy = srodek;
    }

    printf("Pierwiastek nie osiagnieto po %d iteracji, ostatnie przyblizenie to: %.16f \n", limit, srodek);
}


void metoda_newtona(char p, int limit, double tolx, double tolf) {
    double xk = 1.0;

    for (int iter = 0; iter < limit; iter++) {
        double fp = oblicz_f_prim(p, xk);

        if (fp == 0.0) {
            cout << "Pochodna jest bliska zeru, Metoda Newtona jest Rozbieżna" << endl;
            return;
        }

        double xk1   = xk - oblicz_f(p, xk) / fp;
        double blad  = fabs(xk1 - xk);
        double res   = fabs(oblicz_f(p, xk1));

        printf("I: %d xn: %.16f xn+1: %.16f en: %.16f f(xn) %.16f \n",
               iter + 1, xk, xk1, blad, res);

        if (blad <= tolx && res <= tolf) {
            printf("Oczekiwana wartosc uzyskana, konczymy iteracje. Wynik to: %.16f \n", xk1);
            return;
        }

        xk = xk1;
    }

    printf("Pierwiastek nie osiagnieto po %d iteracji, ostatnie przyblizenie to: %.16f", limit, xk);
}


void metoda_siecznych(char p, int limit, double tolx, double tolf) {
    double prev  = -0.5;
    double curr  = 0.5;
    double nxt;

    for (int iter = 0; iter < limit; iter++) {
        double fc   = oblicz_f(p, curr);
        double fp   = oblicz_f(p, prev);
        double mian = (fc - fp) / (curr - prev);

        if (fabs(mian) == 0.0) {
            cout << "Blad: Dzielenie przez zero w metodzie siecznych!" << endl;
            return;
        }

        nxt          = curr - fc / mian;
        double blad  = fabs(nxt - curr);
        double res   = fabs(oblicz_f(p, nxt));

        printf("I: %d xn: %.14f xn1: %.14f xn2: %.14f en: %.14f f(xn) %.14f \n",
               iter + 1, prev, curr, nxt, blad, res);

        if (blad <= tolx && res <= tolf) {
            printf("Oczekiwana wartosc uzyskana, konczymy iteracje. Wynik to: %.16f \n", nxt);
            return;
        }

        prev = curr;
        curr = nxt;
    }

    printf("Pierwiastek nie osiagnieto po %d iteracji, ostatnie przyblizenie to: %.16f", limit, nxt);
}


int main(int argc, char *argv[]) {
    char wybor;
    cout << "Wybierz podpunkt (a lub b): ";
    cin >> wybor;

    if (wybor != 'a' && wybor != 'b') {
        cout << "Podano niepoprawny wybor";
        return 1;
    }

    const int    LIMIT = 50;
    const double TOLX  = 1e-11;
    const double TOLF  = 1e-11;

    cout << "METODA PICARDA" << endl;
    metoda_picarda(wybor, LIMIT, TOLX, TOLF);
    cout << endl;

    cout << "METODA BISEKCJI" << endl;
    metoda_bisekcji(wybor, LIMIT, TOLX, TOLF);
    cout << endl;

    cout << "METODA NEWTONA" << endl;
    metoda_newtona(wybor, LIMIT, TOLX, TOLF);
    cout << endl;

    cout << "METODA SIECZNYCH" << endl;
    metoda_siecznych(wybor, LIMIT, TOLX, TOLF);

    return 0;
}