#include <iostream>
#include <typeinfo>
#include "linalg/Matrix.h"
#include "linalg/Vec.h"

using namespace std;
using namespace geom;

int main_mtxdebug(int argc, char **argv) {
    typedef SimpleMatrix<double,3,4> Ma;
    typedef SimpleMatrix<double,4,2> Mb;
    typedef SimpleMatrix<double,0,0> Mc;
    typedef SimpleMatrix<double,3,2> Md;
    typedef SimpleMatrix<double,3,3> SqM;
    typedef DiagMatrix<double,3,3> diagM;
    typedef PermutationMatrix<3> perM;
    
    cout << IS_MATRIX(Ma) << endl;
    cout << MATRIX_DIM_AGREE(Ma, Ma) << endl;
    cout << MATRIX_MUL_DIM_AGREE(Ma,Mb) << endl;
    cout << MATRIX_MUL_DIM_AGREE(Ma,Ma) << endl;
    cout << MATRIX_MUL_DIM_AGREE(Ma,Mc) << endl;
    cout << typeid(detail::_ImplMtxMul<Ma,Mb>::return_t).name() << endl;
    
    Ma ma;
    Mb mb;
    Md md;
    Mc mc(3,3);
    diagM diag;
    perM per;
    SqM sq;
    
    cout << (ma * mb) << endl;
    
    typedef Vec<double,3> V3;
    typedef Vec<double,4> V4;
    V4 v4(1,2,3,4);
    V3 v3;
    
    cout << MATRIX_MUL_DIM_AGREE(Ma,V4) << endl;
    cout << MATRIX_MUL_DIM_AGREE(Mc,V4) << endl;
    cout << MATRIX_MUL_DIM_AGREE(Ma,V3) << endl;
    
    // static mtx mul dest
    mul(&md, ma, mb);
    cout << md << endl;
    
    // mtx * vec mul
    cout << mul(ma, v4) << endl;
    
    // mtx * vec mul, dest
    mul(&v3, ma, v4);
    cout << v3 << endl;
    
    cout << MATRIX_MUL_DIM_AGREE(diagM, Ma) << endl;
    cout << detail::_ImplMtxAdaptor<SqM, detail::ORIENT_VEC_ROW>::COLDIM << endl;
    
    cout << (diag * ma) << endl;
    
    cout << (diag * v3) << endl;
    
    cout << mul(&v3, diag, v3) << endl;
    
    cout << (per * v3) << endl;
    
    cout << (mc * v3) << endl;
    
    cout << (mc == sq) << endl;
    
    cout << mul(&v3, mc, v3) << endl;
    
    return 0;
}

/* need: vec/diag

class _ImplMtxMul {
class _ImplMtxMul<Mat, DiagMatrix<T,M,N>,                 REQUIRE_MATRIX_T(Mat)> {
class _ImplMtxMul<DiagMatrix<T,M,N>, Mat,                 REQUIRE_MATRIX_T(Mat)> {
class _ImplMtxMul<DiagMatrix<T,L,M1>, DiagMatrix<T,M2,N>, void> {
class _ImplMtxMul<Mat, Vec<T,N>,                          REQUIRE_MATRIX_T(Mat)> {
class _ImplMtxMul<Vec<T,M>, Mat,                          REQUIRE_MATRIX_T(Mat)> { //xx
class _ImplMtxMul<Mat, PermutationMatrix<N>,              REQUIRE_MATRIX_T(Mat)> {
class _ImplMtxMul<PermutationMatrix<M>, Mat,              REQUIRE_MATRIX_T(Mat)> {
class _ImplMtxMul<PermutationMatrix<N1>, Vec<T,N2>,       void> {                  //xx
class _ImplMtxMul<Vec<T,N1>, PermutationMatrix<N2>,       void> {
class _ImplMtxMul<PermutationMatrix<N>, PermutationMatrix<N>, void> {

unit test: 
  - mul every matrix with every other
    - match / mismatch
    - dynamic match / mismatch
    - compare results to SimpleMatrix mul
    - profile speed
  - inv of each type
    - alloc vs. into
    - size 2->6
    - stability / accuracy:
      - divergence of iterated inv
      - inv * orig == ident
      - inv(inv(M)) == orig
    
  - speed tests of other ops

 */
