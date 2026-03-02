/*
 * ============================================================
 *  Laboratorium 2 - Numeryczne obliczanie funkcji
 *  f(x) = x^3 / { 6 * [sinh(x) - x] }   dla x in [10^-10, 10^3]
 * ============================================================
 *
 *  KOMPILACJA:
 *      g++ -O2 -Wall -o lab2 lab2.cpp -lm
 *
 *  URUCHOMIENIE:
 *      ./lab2 dane_do_laboratorium_2.txt
 *
 *  PLIKI WYJSCIOWE:
 *      wyniki.dat   - log10(x), f_naiwny, f_ulepszony, f_ref
 *      bledy.dat    - log10(x), log10(|err_naiwny|), log10(|err_ulepszony|)
 *
 * ============================================================
 *  PROBLEM: KATASTROFALNE ODEJMOWANIE
 * ============================================================
 *
 *  Dla malych x zachodzi:
 *      sinh(x) = x + x^3/6 + x^5/120 + x^7/5040 + ...
 *
 *  Wiec:  sinh(x) - x  ~  x^3/6  (bardzo mala liczba gdy x -> 0)
 *
 *  Przyklad: x = 0.001
 *    sinh(x) = 0.001 000 000 166 666 750 ...   (16 cyfr znaczacych)
 *          x = 0.001 000 000 000 000 000 ...
 *  roznica  = 0.000 000 000 166 666 750 ...   <- tylko 10 cyfr zostalo!
 *                              ^^^^^^^^^
 *                        6 cyfr bezpowrotnie utraconych
 *
 *  Im mniejsze x, tym wiecej cyfr ginie.
 *  Dla x < 1e-8:  sinh(x) - x zaokragla sie do DOKLADNIE 0
 *                 -> dzielenie przez zero -> wynik = +Inf
 *
 * ============================================================
 *  ROZWIAZANIE: SZEREG TAYLORA
 * ============================================================
 *
 *  Rozwijamy mianownik:
 *    sinh(x) - x = x^3/3! + x^5/5! + x^7/7! + x^9/9! + ...
 *                = x^3/6  + x^5/120 + x^7/5040 + x^9/362880 + ...
 *
 *  Wyciagamy x^3/6 przed nawias:
 *    sinh(x) - x = (x^3/6) * [1 + x^2/20 + x^4/840 + x^6/60480 + x^8/6652800 + ...]
 *
 *  Wspolczynniki (skad sie biora?):
 *    x^5/120  / (x^3/6)  = 6*x^2/120  = x^2/20
 *    x^7/5040 / (x^3/6)  = 6*x^4/5040 = x^4/840
 *    x^9/362880/(x^3/6)  = 6*x^6/362880= x^6/60480
 *    x^11/39916800/(x^3/6)= 6*x^8/39916800 = x^8/6652800
 *
 *  Podstawiamy do f(x):
 *    f(x) = x^3 / { 6 * (x^3/6) * [1 + x^2/20 + ...] }
 *         = 1   / [1 + x^2/20 + x^4/840 + x^6/60480 + x^8/6652800 + ...]
 *
 *  Brak odejmowania bliskich liczb -> blad maszynowy wszedzie ~2e-16 !
 *
 *  PROG PRZELACZANIA: |x| < 0.1
 *    Dla x = 0.1: nastepny wyraz szeregu ~ x^10/... ~ 1e-17 < eps_mach
 *    Dla x >= 0.1: wzor naiwny stabilny (sinh(x) >> x, nie ma bliskiego odejmowania)
 *
 *  TYP ZMIENNYCH: double (64-bit IEEE 754)
 *    - 53 bity mantysy  ~ 15-16 cyfr dziesietnych
 *    - blad maszynowy eps_mach ~ 2.2e-16
 *    - float (32-bit, ~7 cyfr) - za malo precyzyjny
 *    - long double (80-bit, ~19 cyfr) - zbedny przy poprawnym algorytmie
 * ============================================================
 */

#include <cstdio>
#include <cmath>
#include <cstdlib>

/* ----------------------------------------------------------
 *  Algorytm 1: NAIWNY
 *  Bezposredni wzor ze standardowa funkcja sinh()
 *  UWAGA: duzy blad dla malych x !
 * ---------------------------------------------------------- */
double f_naiwny(double x)
{
    return (x * x * x) / (6.0 * (sinh(x) - x));
}

/* ----------------------------------------------------------
 *  Algorytm 2: ULEPSZONY
 *  - dla |x| < 0.1: szereg Taylora (bez odejmowania)
 *  - dla |x| >= 0.1: wzor naiwny (numerycznie stabilny)
 * ---------------------------------------------------------- */
double f_ulepszony(double x)
{
    if (fabs(x) < 0.1)
    {
        /*
         * f(x) = 1 / (1 + x^2/20 + x^4/840 + x^6/60480 + x^8/6652800)
         *
         * Obliczanie schematem Hornera (minimalna liczba mnozen):
         * mianownik = 1 + x2*(A + x2*(B + x2*(C + x2*D)))
         */
        double x2 = x * x;
        double mian = 1.0 + x2 * ( 1.0/20.0
                          + x2 * ( 1.0/840.0
                          + x2 * ( 1.0/60480.0
                          + x2 *   1.0/6652800.0 )));
        return 1.0 / mian;
    }
    else
    {
        /* Wzor bezposredni - stabilny dla wiekszych x */
        return f_naiwny(x);
    }
}

/* ----------------------------------------------------------
 *  Pomocnicza: log10 z bezwzglednego bledu wzglednego
 *  Zwraca -99 gdy blad jest dokladnie 0
 * ---------------------------------------------------------- */
static double log10_err(double obliczone, double referencyjne)
{
    if (referencyjne == 0.0) return -99.0;
    double err = fabs((obliczone - referencyjne) / referencyjne);
    if (err == 0.0) return -99.0;
    return log10(err);
}

/* ========================================================== */
int main(int argc, char *argv[])
{
    const char *nazwa = (argc >= 2) ? argv[1] : "dane_do_laboratorium_2.txt";

    FILE *fin = fopen(nazwa, "r");
    if (!fin) {
        fprintf(stderr, "BLAD: Nie mozna otworzyc '%s'\n", nazwa);
        return 1;
    }

    FILE *fw = fopen("wyniki.dat", "w");
    FILE *fb = fopen("bledy.dat",  "w");
    if (!fw || !fb) {
        fprintf(stderr, "BLAD: Nie mozna otworzyc plikow wyjsciowych.\n");
        return 1;
    }

    fprintf(fw, "# %12s  %20s  %20s  %20s\n",
            "log10(x)", "f_naiwny", "f_ulepszony", "f_referencyjne");
    fprintf(fb, "# %12s  %20s  %20s\n",
            "log10(x)", "log10|err_naiwny|", "log10|err_ulepszony|");

    /* Pomin 3 linie naglowka pliku danych */
    char linia[512];
    for (int i = 0; i < 3; i++) {
        if (!fgets(linia, sizeof(linia), fin)) break;
    }

    double log10x, x, f_ref;
    int n_inf = 0, n_ok = 0;

    while (fscanf(fin, "%lf %lf %lf", &log10x, &x, &f_ref) == 3)
    {
        double fn = f_naiwny(x);
        double fi = f_ulepszony(x);

        fprintf(fw, "  %12.5f  %20.10e  %20.10e  %20.10e\n",
                log10x, fn, fi, f_ref);

        double lerr_n = log10_err(fn, f_ref);
        double lerr_i = log10_err(fi, f_ref);

        /* Zapisz bledy tylko gdy obie wartosci sa skonczONE */
        if (std::isfinite(lerr_n) && std::isfinite(lerr_i)) {
            fprintf(fb, "  %12.5f  %20.6f  %20.6f\n",
                    log10x, lerr_n, lerr_i);
            n_ok++;
        } else {
            /* Naiwny daje Inf - zaznacz to jako blad = 0 (wynik nieskoncz.) */
            n_inf++;
        }
    }

    fclose(fin); fclose(fw); fclose(fb);

    printf("\n");
    printf("  Plik wejsciowy : %s\n",  nazwa);
    printf("  Punkty ok      : %d\n",  n_ok);
    printf("  Punkty Inf/NaN : %d  (naiwny: dzielenie przez ~0)\n", n_inf);
    printf("  --> wyniki.dat, bledy.dat\n");
    printf("  Wykres: gnuplot wykres.gnu\n\n");
    return 0;
}