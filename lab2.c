#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PROG 0.1L

long double fx_naiwnie(long double x){
    long double x3 = x*x*x;
    long double sinh = sinhl(x);
    long double roznica = sinh - x;
    long double mianownik = 6.0L * roznica;
    return x3 / mianownik;
}

long double fx_taylor(long double x) {
    long double x2 = x * x;
    long double x4 = x2 * x2;
    long double x6 = x4 * x2;
    long double x8 = x6 * x2;

    long double skladnik1 = x2 / 20.0L;
    long double skladnik2 = x4 / 840.0L;
    long double skladnik3 = x6 / 60480.0L;
    long double skladnik4 = x8 / 6652800.0L;

    long double mianownik = 1.0L + skladnik1 + skladnik2 + skladnik3 + skladnik4;

    return 1.0L / mianownik;
}

long double blad(long double dokl, long double obl) {
    return fabsl(obl - dokl)/ fabsl(dokl);
}

int licz_wiersze(const char* nazwa) {
    FILE *f = fopen(nazwa, "r");
    if (!f) return -1;

    char bufor[256];
    for (int i = 0;i<3;i++) fgets(bufor, sizeof(bufor),f);

    int n = 0;
    long double a,b,c;
    while (fscanf(f,"%Lg" "%Lg" "%Lg", &a,&b,&c) == 3) n++;
    fclose(f);
    return n;
}

int wczytaj(const char* nazwa, long double* log10x, long double* x, long double* fx_dokl) {
    FILE *f = fopen(nazwa, "r");
    if (!f) return -1;
    char bufor[256];
    for (int i = 0;i<3;i++) fgets(bufor, sizeof(bufor),f);
    int n = 0;
    while (fscanf(f,"%Lg" "%Lg" "%Lg", &log10x[n], &x[n], &fx_dokl[n]) == 3) n++;
   fclose(f);
    return n;
}

int main(void) {
    const char* plik_danych = "dane_do_laboratorium_2.txt";

    int n = licz_wiersze(plik_danych);
    if (n <= 0) {
        fprintf(stderr, "brak pliku danych\n");
        return -1;
    }
    long double* log10x = malloc(n*sizeof(long double));
    long double* x = malloc(n*sizeof(long double));
    long double* fx_dokl = malloc(n*sizeof(long double));

    wczytaj(plik_danych, log10x, x, fx_dokl);
    printf("wczytano %d punktow\n", n);

    FILE* out1 = fopen("bledy_podstawowa.txt", "w");
    fprintf(out1, "# log10x  log10_blad\n");
    for (int i = 0; i < n; i++) {
        long double b = blad(fx_dokl[i], fx_naiwnie(x[i]));
        fprintf(out1, "%.6Lg  %.6Lg\n", log10x[i], b > 0.0L ? log10l(b) : -20.0L);
    }
    fclose(out1);

    /* --- metoda ulepszona --- */
    FILE* out2 = fopen("bledy_ulepszona.txt", "w");
    fprintf(out2, "# log10x  log10_blad\n");
    for (int i = 0; i < n; i++) {
        long double fx = (fabsl(x[i]) < PROG) ? fx_taylor(x[i]) : fx_naiwnie(x[i]);
        long double b  = blad(fx_dokl[i], fx);
        fprintf(out2, "%.6Lg  %.6Lg\n", log10x[i], b > 0.0L ? log10l(b) : -20.0L);
    }
    fclose(out2);

    printf("Zapisano: bledy_podstawowa.txt, bledy_ulepszona.txt\n");

    free(log10x);
    free(x);
    free(fx_dokl);
    return 0;
}
