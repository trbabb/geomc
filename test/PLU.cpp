#include <iostream>
#include "random/MTRand.h"
#include "random/RandomTools.h"
#include "linalg/Matrix.h"

using namespace std;
using namespace geom;

int poo(int argc, char **argv) {
//int main(int argc, char **argv) {
    const static index_t sz = 5;
    SimpleMatrix<double,sz,sz> m;
    
    Random *rng = new MTRand(20928234);
    
    for (index_t r = 0; r < m.rows(); r++) {
        for (index_t c = 0; c < m.cols(); c++) {
            m.set(r,c,rng->rand<double>());
        }
    }
    
    PLUDecomposition<double,sz,sz> plu = plu_decompose(m);
    
    cout << "singular? " << plu.isSingular() << endl;
    
    cout << m << endl;
    cout << "P = " << endl << plu.P << endl;
    cout << "L = " << endl << plu.L << endl;
    cout << "U = " << endl << plu.U << endl;
    
    cout << (plu.L * plu.U) << endl;
    cout << (plu.P * m) << endl;
    
    Vec<double,sz> b;
    Vec<double,sz> x;
    for (index_t i = 0; i < sz; i++) {
        b[i] = rng->rand<double>();
    }
    
    plu.linearSolve<double>(x.begin(), b.begin());
    
    cout << "b  = " << b << endl;
    cout << "x  = " << x << endl;
    cout << "Mx = " << (m * x) << endl;
    cout << "residual: " << (m * x) - b << endl;
    
    SimpleMatrix<double,sz,sz> m_inv;
    plu.inverse(&m_inv);
    
    bool ok = true;
    
    cout << "inv: " << endl;
    cout << m_inv << endl;
    cout << "M * M^-1:" << endl << (m * m_inv) << endl;
    cout << "inv(M): " << endl << m * inv(m, &ok) << endl;
    
    return 0;
}

/*
PA = LU
Pb = LUx
Pb = Ly
 y = Ux

*/
