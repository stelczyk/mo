#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>

#define M_PI 3.14159265358979323846

template<typename T>
T func(T x) { return std::cos(x); }

template<typename T>
T exact_deriv(T x) { return -std::sin(x); }


template<typename T>
T forward2(T x, T h) { return (func(x + h) - func(x)) / h; }

template<typename T>
T forward3(T x, T h) {
    return (-T(3)/T(2)*func(x) + T(2)*func(x+h) - T(1)/T(2)*func(x+2*h)) / h;
}

template<typename T>
T central2(T x, T h) { return (func(x + h) - func(x - h)) / (T(2)*h); }

template<typename T>
T backward2(T x, T h) { return (func(x) - func(x - h)) / h; }

template<typename T>
T backward3(T x, T h) {
    return (T(3)/T(2)*func(x) - T(2)*func(x-h) + T(1)/T(2)*func(x-2*h)) / h;
}

template<typename T>
void compute_errors(
    const std::string& label,
    T x_point,
    T (*approx)(T, T),
    std::ofstream& out,
    int n_steps = 200
) {
    T exact = exact_deriv(x_point);
    for (int i = 0; i <= n_steps; ++i) {
        double log_h = -20.0 + 19.0 * i / n_steps;
        T h = static_cast<T>(std::pow(10.0, log_h));
        T approx_val = approx(x_point, h);
        T err = std::abs(approx_val - exact);

        if (err > T(0)) {
            out << log_h << " " << std::log10(static_cast<double>(err)) << " " << label << "\n";
        }
    }
    out << "\n\n";
}

template<typename T>
void print_order(const std::string& name, T x_point, T (*approx)(T, T), int theoretical_p) {
    // Dobieramy kroki h, dla których dominuje błąd metody (błąd obcięcia),
    // a nie błąd reprezentacji zmiennoprzecinkowej (maszynowy).
    T h1 = 1e-1;
    T h2 = 1e-3;

    T exact = exact_deriv(x_point);
    T err1 = std::abs(approx(x_point, h1) - exact);
    T err2 = std::abs(approx(x_point, h2) - exact);

    double p_exp = 0.0;
    if (err1 > 0 && err2 > 0) {
        p_exp = (std::log10(static_cast<double>(err1)) - std::log10(static_cast<double>(err2))) /
                (std::log10(static_cast<double>(h1)) - std::log10(static_cast<double>(h2)));
    }

    std::cout << std::left << std::setw(25) << name
              << " | Teoretyczny p: " << theoretical_p
              << " | Eksperymentalny p: " << std::fixed << std::setprecision(4) << p_exp << "\n";
}


void generate_and_run_gnuplot() {
    std::ofstream script("plot_errors.plt");

    script << "set terminal pngcairo size 1100,750 enhanced font 'Arial,11'\n";
    script << "set grid\n";
    script << "set xlabel 'log_{10}(h)'\n";
    script << "set ylabel 'log_{10}(|blad bezwzgledny|)'\n";
    script << "set key outside right top\n\n";

    // Definicja stylów: Double (linia ciągła), Long Double (linia przerywana 'dt 2')
    script << "set style line 1 lc rgb 'blue' lw 2\n";    // Metoda 1 (np. fwd2 / bwd2)
    script << "set style line 2 lc rgb 'red' lw 2\n";      // Metoda 2 (np. fwd3 / bwd3)
    script << "set style line 3 lc rgb 'dark-green' lw 2\n"; // Centralna

    // 1. Wykres dla x0 = 0
    script << "set output 'wykres_x0_porownanie.png'\n";
    script << "set title 'Porownanie Double vs Long Double: x = 0'\n";
    script << "plot 'errors_double.dat' index 0 with lines ls 1 title 'Dbl: fwd 2-pkt', \\\n";
    script << "     'errors_double.dat' index 1 with lines ls 2 title 'Dbl: fwd 3-pkt', \\\n";
    script << "     'errors_longdouble.dat' index 0 with lines ls 1 dt 2 title 'LDbl: fwd 2-pkt', \\\n";
    script << "     'errors_longdouble.dat' index 1 with lines ls 2 dt 2 title 'LDbl: fwd 3-pkt'\n\n";

    // 2. Wykres dla x_mid = pi/4
    script << "set output 'wykres_x_mid_porownanie.png'\n";
    script << "set title 'Porownanie Double vs Long Double: x = pi/4'\n";
    script << "plot 'errors_double.dat' index 2 with lines ls 3 title 'Dbl: central 2-pkt', \\\n";
    script << "     'errors_longdouble.dat' index 2 with lines ls 3 dt 2 title 'LDbl: central 2-pkt'\n\n";

    // 3. Wykres dla xn = pi/2
    script << "set output 'wykres_xn_porownanie.png'\n";
    script << "set title 'Porownanie Double vs Long Double: x = pi/2'\n";
    script << "plot 'errors_double.dat' index 3 with lines ls 1 title 'Dbl: bwd 2-pkt', \\\n";
    script << "     'errors_double.dat' index 4 with lines ls 2 title 'Dbl: bwd 3-pkt', \\\n";
    script << "     'errors_longdouble.dat' index 3 with lines ls 1 dt 2 title 'LDbl: bwd 2-pkt', \\\n";
    script << "     'errors_longdouble.dat' index 4 with lines ls 2 dt 2 title 'LDbl: bwd 3-pkt'\n";

    script.close();
    std::system("gnuplot plot_errors.plt");
}

int main() {
    const double x_left  = 0.0;
    const double x_mid   = M_PI / 4.0;
    const double x_right = M_PI / 2.0;

    // --- Wypisywanie rzędów dokładności ---
    std::cout << "--- RZEDY DOKLADNOSCI (DOUBLE) ---\n";
    print_order<double>("x0_fwd2 (2-punktowe)", x_left,  forward2<double>, 1);
    print_order<double>("x0_fwd3 (3-punktowe)", x_left,  forward3<double>, 2);
    print_order<double>("xm_cen2 (centralne)",  x_mid,   central2<double>, 2);
    print_order<double>("xn_bwd2 (2-punktowe)", x_right, backward2<double>, 1);
    print_order<double>("xn_bwd3 (3-punktowe)", x_right, backward3<double>, 2);

    std::cout << "\n--- RZEDY DOKLADNOSCI (LONG DOUBLE) ---\n";
    print_order<long double>("x0_fwd2 (2-punktowe)", (long double)x_left,  forward2<long double>, 1);
    print_order<long double>("x0_fwd3 (3-punktowe)", (long double)x_left,  forward3<long double>, 2);
    print_order<long double>("xm_cen2 (centralne)",  (long double)x_mid,   central2<long double>, 2);
    print_order<long double>("xn_bwd2 (2-punktowe)", (long double)x_right, backward2<long double>, 1);
    print_order<long double>("xn_bwd3 (3-punktowe)", (long double)x_right, backward3<long double>, 2);
    std::cout << "\n";

    // --- Generowanie plików z danymi ---
    // Generowanie errors_double.dat
    {
        std::ofstream out("errors_double.dat");
        out << std::scientific << std::setprecision(15);
        compute_errors<double>("x0_fwd2", x_left,  forward2<double>,  out);
        compute_errors<double>("x0_fwd3", x_left,  forward3<double>,  out);
        compute_errors<double>("xm_cen2", x_mid,   central2<double>,  out);
        compute_errors<double>("xn_bwd2", x_right, backward2<double>, out);
        compute_errors<double>("xn_bwd3", x_right, backward3<double>, out);
    }

    // Generowanie errors_longdouble.dat
    {
        std::ofstream out("errors_longdouble.dat");
        out << std::scientific << std::setprecision(18);
        compute_errors<long double>("x0_fwd2", (long double)x_left,  forward2<long double>,  out);
        compute_errors<long double>("x0_fwd3", (long double)x_left,  forward3<long double>,  out);
        compute_errors<long double>("xm_cen2", (long double)x_mid,   central2<long double>,  out);
        compute_errors<long double>("xn_bwd2", (long double)x_right, backward2<long double>, out);
        compute_errors<long double>("xn_bwd3", (long double)x_right, backward3<long double>, out);
    }

    generate_and_run_gnuplot();
    return 0;
}