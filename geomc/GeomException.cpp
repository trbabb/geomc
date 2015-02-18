/*
 * GeometryException.cpp
 *
 *  Created on: Feb 22, 2009
 *      Author: Tim Babb
 */

#include <stdexcept>
#include <stdio.h>
#include <geomc/GeomException.h>

/////////////// GeomException ///////////////

GeomException::GeomException(const std::string& msg) : std::runtime_error(msg) {
    //do nothing
}

GeomException::~GeomException() throw() {
    //do nothing
}

/////////////// DimensionMismatch ///////////////

DimensionMismatchException::DimensionMismatchException(index_t a_0, index_t a_1, index_t b_0, index_t b_1) throw () :
        std::runtime_error("dimension mismatch"),
        a_0(a_0),
        a_1(a_1),
        b_0(b_0),
        b_1(b_1) { /* do nothing */ }

DimensionMismatchException::~DimensionMismatchException() throw() { /* do nothing */ }

/////////////// NonsquareMatrix ///////////////

NonsquareMatrixException::NonsquareMatrixException(index_t rows, index_t cols) throw () :
        std::runtime_error("nonsquare matrix"),
        rows(rows), cols(cols) { /* do nothing */ }

NonsquareMatrixException::~NonsquareMatrixException() throw() { /* do nothing */ }