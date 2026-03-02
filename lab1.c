

#include <stdio.h>
#include <math.h>

void analyze_float() {
    float eps = 1.0f;
    int bits = 0;


    while ((float)(1.0f + (eps / 2.0f)) > 1.0f) {
        eps /= 2.0f;
        bits++;
    }

    printf("FLOAT\n");
    printf("Epsilon maszynowy: %.20e\n", eps);
    printf("Liczba bitow mantysy: %d (wliczajac bit ukryty: %d)\n", bits, bits + 1);
    printf("Cyfry znaczace (log10(1/eps)): ~%.2f\n\n", log10(1.0/eps));
}

void analyze_double() {
    double eps = 1.0;
    int bits = 0;

    while ((double)(1.0 + (eps / 2.0)) > 1.0) {
        eps /= 2.0;
        bits++;
    }

    printf("DOUBLE\n");
    printf("Epsilon maszynowy: %.20e\n", eps);
    printf("Liczba bitow mantysy: %d (wliczajac bit ukryty: %d)\n", bits, bits + 1);
    printf("Cyfry znaczace (log10(1/eps)): ~%.2f\n\n", log10(1.0/eps));
}

void analyze_long_double() {
    long double eps = 1.0L;
    int bits = 0;

    while ((long double)(1.0L + (eps / 2.0L)) >  1.0L) {
        eps /= 2.0L;
        bits++;
    }

    printf("LONG DOUBLE\n");
    printf("Epsilon maszynowy: %.25Le\n", eps);
    printf("Liczba bitow mantysy: %d (wliczajac bit ukryty: %d)\n", bits, bits + 1);
    printf("Cyfry znaczace: ~%.2f\n\n", (double)log10l(1.0L/eps));
}

int main() {
    analyze_float();
    analyze_double();
    analyze_long_double();

    return 0;
}