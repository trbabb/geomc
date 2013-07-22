/*
 * Matrix.h
 * 
 * There are currently six distinct matrix template classes.
 * They are all interoperable, and provide iterators which
 * are functionally interchangeable with each other and
 * with pointers (row iterators, column iterators, region
 * iterators, and row-major matrix body iterators).
 * 
 * These iterators are compatible with std::copy(). Some
 * of them are not writeable (or may possibly throw an
 * error upon writing) in the case where a matrix element
 * does not have a corresponding memory location (DiagonalMatrix
 * is one such case, for example; off-diagonals are not stored).
 * 
 * In general, this scheme was chosen to satisfy the following
 * requirements:
 * 
 * - All matrix types must interoperate via mul() and copy().
 * - Element access shall be non-virtual and inline-able wherever possible.
 * - Copies to and from contiguous memory shall be fast where the
 *   internal matrix representation is contiguous.
 * - Handles to matrices of arbitrary type are possible.
 * 
 *  Created on: Oct 27, 2010
 *      Author: tbabb
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include "LinalgTypes.h"

// matrix types
// each includes matrixdetail and matrixglue

#include "mtxtypes/SimpleMatrix.h"
#include "mtxtypes/AugmentedMatrix.h"
#include "mtxtypes/DiagMatrix.h"
#include "mtxtypes/SparseMatrix.h"
#include "mtxtypes/GeneralMatrix.h"
#include "mtxtypes/PermutationMatrix.h"

// includes matrixfunctionimpl
#include "mtxdetail/MatrixFunction.h"

//TODO: templatize memory layout choice
//TODO: reduce bloat, particularly in matrix inv case.
//TODO: determinant
//TODO: matrix kernel (ND normal-finder)
//TODO: clean arbitrary Matrix handle construction
//TODO: verify correct matrix template function resolution.

//future: function- or expression-valued matrices?
//future: look for ways to keep object code bloat down.
//future: mechanism for col-major mtxs? (iterator based wrapper?) 
//        (could clutter code, because then other classes must be similarly templated).

//refactoring:

//TODO: copy methods for subsets that can grab chunks by row if mem is contiguous.
//TODO: iterators for nonzero entries.
//TODO: set() should return reference to (or value of) new element, so that `z = (mtx[x][y] = foo)` can
//      evaluate properly with proxy references.

#endif /* MATRIX_H_ */
