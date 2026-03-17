#include <stdio.h>

void analizuj_float() {
    float eps = 1.0f;
    int bity = 0;

    while (1.0f + (eps / 2.0f) > 1.0f) {
        eps = eps / 2.0f;
        bity = bity + 1;
    }



    printf("--- Typ FLOAT ---\n");
    printf("Liczba bitow mantysy: %d\n", bity);
    printf("Cyfry znaczace: 15\n");
    printf("Epsylon: %.100e\n\n", eps);
}

void analizuj_double() {
    double eps = 1.0;
    int bity = 0;

    while (1.0 + (eps / 2.0) > 1.0) {
        eps = eps / 2.0;
        bity = bity + 1;
    }



    printf("--- Typ DOUBLE ---\n");
    printf("Liczba bitow mantysy: %d\n", bity);
    printf("Cyfry znaczace: 25\n");
    printf("Epsylon: %.100e\n\n", eps);
}

void analizuj_long_double() {
    long double eps = 1.0L;
    int bity = 0;

    while (1.0L + (eps / 2.0L) > 1.0L) {
        eps = eps / 2.0L;
        bity = bity + 1;
    }


    printf("--- Typ LONG DOUBLE ---\n");
    printf("Liczba bitow mantysy: %d\n", bity);
    printf("Cyfry znaczace: 44\n");
    printf("Epsylon: %.100Le\n\n", eps);
}

int main() {
    analizuj_float();
    analizuj_double();
    analizuj_long_double();

    return 0;
}