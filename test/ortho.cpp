#include <iostream>
#include <geomc/linalg/Orthogonal.h>
#include <geomc/random/RandomTools.h>

#define N 6

using namespace std;
using namespace geom;

void ortho() {
    Sampler<double> smp;
    Vec<double,N> basis[N-1];
    for (int i = 0; i < N-1; i++) {
        basis[i] = smp.unit<N>();
        cout << basis[i] << "\n";
    }
    cout << "\n";
    Vec<double,N> o = orthogonal(basis);
    cout << "N: " << o << "\n";
    for (int i = 0; i < N-1; i++) {
        cout << "dot: " << basis[i].dot(o) << "\n";
    }
}


void nullsp() {
    Sampler<double> smp;
    Vec<double,N> basis[N];
    index_t n = getRandom()->rand<index_t>(N-2) + 1;
    cout << "given basis:\n";
    for (int i = 0; i < n; i++) {
        basis[i] = smp.unit<N>();
        cout << basis[i] << "\n";
    }
    
    nullspace(basis, n, basis + n);
    
    cout << "\nnullspace basis:\n";
    for(int i = 0; i < N-n; i++) {
        cout << basis[n+i] << "\n";
    }
    cout << "\ndots with original basis:\n";
    for (int b = 0; b < N-n; b++) {
        for (int q = 0; q < n; q++) {
            cout << "dot: " << basis[n+b].dot(basis[q]) << "\n";
        }
    }
    cout << "\ndots with selves:\n";
    for (int b = 0; b < N-n; b++) {
        for (int b2 = 0; b2 < N-n; b2++) {
            if (b2 == b) continue;
            cout << "dot: " << basis[n+b].dot(basis[n+b2]) << "\n";
        }
    }
    cout << "\nlengths:\n";
    for(int b = 0; b < N-n; b++) {
        cout << basis[n+b].mag() << "\n";
    }
}

int main(int argc, char **argv) {
    cout << "\n====== orthogonal ======\n";
    ortho();
    cout << "\n\n====== nullspace ======\n";
    nullsp();
    return 0;
}
