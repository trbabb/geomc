/*
 * GeometryException.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#ifndef GEOMETRYEXCEPTION_H_
#define GEOMETRYEXCEPTION_H_

#include <geomc/geomc_defs.h>
 
namespace geom {

///////////////////////

class GeomException {
public:
    GeomException(const char* msg) throw () : msg(msg) {}

    const char* what() { return msg; }

    const char* msg;
};

///////////////////////

//todo: this should be factored to linalg, and use MatrixDim
//todo: rename MatrixDim to Size2d? make a Size<#>::type?

class DimensionMismatchException : public GeomException {
public:
    index_t a_0, a_1, b_0, b_1;

    DimensionMismatchException(index_t a_0, index_t a_1, index_t b_0, index_t b_1) throw () :
        GeomException("dimension mismatch"),
        a_0(a_0),
        a_1(a_1),
        b_0(b_0),
        b_1(b_1) {}
};

///////////////////////

class NonsquareMatrixException : public GeomException {
public:
    index_t rows, cols;

    NonsquareMatrixException(index_t rows, index_t cols) throw () :
        GeomException("nonsquare matrix"),
        rows(rows), cols(cols) {}
};

} // namespace geom

#endif /* GEOMETRYEXCEPTION_H_ */
