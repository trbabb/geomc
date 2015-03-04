/*
 * GeometryException.h
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#ifndef GEOMETRYEXCEPTION_H_
#define GEOMETRYEXCEPTION_H_

#include <stdexcept>
#include <geomc/geomc_defs.h>
 
namespace geom {

///////////////////////

class GeomException : virtual public std::runtime_error {
public:
    GeomException(const std::string& msg);
    virtual ~GeomException() throw ();
};

///////////////////////

//todo: this should be factored to linalg, and use MatrixDim
//todo: rename MatrixDim Coords2d

class DimensionMismatchException : virtual public std::runtime_error {
public:
    index_t a_0, a_1, b_0, b_1;
    
    DimensionMismatchException(index_t a_0, index_t a_1, index_t b_0, index_t b_1) throw ();
    virtual ~DimensionMismatchException() throw ();
};

///////////////////////

class NonsquareMatrixException : virtual public std::runtime_error {
public:
    index_t rows, cols;
    
    NonsquareMatrixException(index_t rows, index_t cols) throw ();
    virtual ~NonsquareMatrixException() throw ();
};

} // namespace geom

#endif /* GEOMETRYEXCEPTION_H_ */
