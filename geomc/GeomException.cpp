/*
 * GeometryException.cpp
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#include <geomc/GeomException.h>


namespace geom {


/////////////// GeomException ///////////////

GeomException::GeomException(const char* msg) throw () : msg(msg) { /* do nothing */ }

const char* GeomException::what() { return msg; }

/////////////// DimensionMismatch ///////////////

DimensionMismatchException::DimensionMismatchException(index_t a_0, index_t a_1, index_t b_0, index_t b_1) throw () :
        GeomException("dimension mismatch"),
        a_0(a_0),
        a_1(a_1),
        b_0(b_0),
        b_1(b_1) { /* do nothing */ }

/////////////// NonsquareMatrix ///////////////

NonsquareMatrixException::NonsquareMatrixException(index_t rows, index_t cols) throw () :
        GeomException("nonsquare matrix"),
        rows(rows), cols(cols) { /* do nothing */ }


} // namespace geom
