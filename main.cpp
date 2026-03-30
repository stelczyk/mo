#include <iostream>
#include <cmath>

using namespace std;

void metoda_newtona_uogolniona() {

    const int    LIMIT = 20;
    const double TOLX  = 1e-11;
    const double TOLF  = 1e-11;

    double x = 1.0, y = 1.0, z = 1.0;

    for (int iter = 0; iter < LIMIT; iter++) {

        double f1 = x*x + y*y + z*z - 4.0;
        double f2 = x*x + y*y / 2.0  - 1.0;
        double f3 = x*y - 0.5;

        // wyznacznik podmacierzy 2x2 powstalej po eliminacji z
        double det = 2.0*x*x - y*y;

        if (fabs(det) == 0.0 || fabs(y) == 0.0 || fabs(z) == 0.0) {
            cout << "Dzielenie przez zero w iteracji " << iter + 1
                 << ". Metoda rozbiezna." << endl;
            return;
        }

        double dy  = (2.0*x*f3 - f2*y) / det;
        double dx = (f3 - x * dy) / y;
        double dz = (f1 - 2.0*x*dx - 2.0*y*dy) / (2.0 * z);

        double est_x = fabs(dx);
        double est_y = fabs(dy);
        double est_z = fabs(dz);
        double res_1 = fabs(f1);
        double res_2 = fabs(f2);
        double res_3 = fabs(f3);
        double est_max = fmax(fmax(est_x, est_y), est_z);
        double res_max = fmax(fmax(res_1, res_2), res_3);

        printf("I: %d\n"
               "   x:  %.14f   y:  %.14f   z:  %.14f\n"
               "   dx: %.14f   dy: %.14f   dz: %.14f\n"
               "   f1: %.14f   f2: %.14f   f3: %.14f\n"
               "   estymator_max: %.14f   residuum_max: %.14f\n",
               iter + 1,
               x, y, z,
               est_x, est_y, est_z,
               res_1, res_2, res_3,
               est_max, res_max);

        bool kryterium_kroku    = est_x <= TOLX && est_y <= TOLX && est_z <= TOLX;
        bool kryterium_residuum = res_1  <= TOLF && res_2  <= TOLF && res_3  <= TOLF;

        if (kryterium_kroku && kryterium_residuum) {
            printf("Osiagnieto rozwiazanie w %d iteracjach.\n"
                   "Wektor wynikowy:\n"
                   "  x = %.17f\n"
                   "  y = %.17f\n"
                   "  z = %.17f\n",
                   iter + 1, x - dx, y - dy, z - dz);
            return;
        }

        x -= dx;
        y -= dy;
        z -= dz;
    }

    printf("Nie osiagnieto zbieznosci po %d iteracjach.\n"
           "Ostatnie przyblizenie:\n"
           "  x = %.17f\n  y = %.17f\n  z = %.17f\n",
           LIMIT, x, y, z);
}

int main() {
    cout << "METODA NEWTONA UOGOLNIONA" << endl;
    cout << "Uklad: x2+y2+z2=4 | x2+y2/2=1 | xy=1/2" << endl;
    cout << "Punkt startowy: x=1, y=1, z=1" << endl << endl;
    metoda_newtona_uogolniona();
    return 0;
}