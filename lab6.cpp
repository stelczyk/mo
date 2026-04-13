#include <iostream>
#include <vector>
using namespace std;

void thomas_dekompozycja(vector<double>& d, const vector<double>& l, const vector<double>& u,
                         vector<double>& gamma) {
    int N = d.size();
    gamma.resize(N);
    gamma[0] = d[0];
    for (int i = 1; i < N; i++)
        gamma[i] = d[i] - (l[i-1] / gamma[i-1]) * u[i-1];
}

void thomas_rozwiaz(const vector<double>& gamma, const vector<double>& l,
                    const vector<double>& u, const vector<double>& b,
                    vector<double>& x) {
    int N = gamma.size();
    vector<double> tau(N);

    tau[0] = b[0];
    for (int i = 1; i < N; i++)
        tau[i] = b[i] - (l[i-1] / gamma[i-1]) * tau[i-1];

    x.resize(N);
    x[N-1] = tau[N-1] / gamma[N-1];
    for (int i = N-2; i >= 0; i--)
        x[i] = (tau[i] - u[i] * x[i+1]) / gamma[i];
}

int main() {

    vector<double> d = {100, 200, 300, 200, 100};
    vector<double> l = {2, 4, -6, -8};
    vector<double> u = {-1, -3, 5, -7};
    vector<double> b = {199, 195, 929, 954, 360};

    vector<double> gamma;
    vector<double> x;

    thomas_dekompozycja(d, l, u, gamma);
    thomas_rozwiaz(gamma, l, u, b, x);

    cout << "Rozwiazanie x:" << endl;
    for (int i = 0; i < (int)x.size(); i++)
        cout << "x[" << i+1 << "] = " << x[i] << endl;

    return 0;
}