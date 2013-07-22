/*
 * MatrixMacro.h
 *
 *  Created on: Jun 23, 2013
 *      Author: tbabb
 */

#ifndef MATRIXMACRO_H_
#define MATRIXMACRO_H_

// TODO: some of these could be made into friendlier templates
//       (i.e. dim agreement), which could clean up template error reporting.
//       this would also return valid results if T is not a matrix, 
//       though beware of places where you are currently counting on 
//       substitution failure.

//////////////////// Matrix template macros ////////////////////

// Mt is a type
#define REQUIRE_MATRIX_T(Mt) typename boost::enable_if<boost::is_base_of<geom::detail::MatrixBase<typename Mt::elem_t, Mt::ROWDIM, Mt::COLDIM, Mt>, Mt> >::type

// Mt is a type
#define IS_MATRIX(Mt) (boost::is_base_of<geom::detail::MatrixBase<typename Mt::elem_t, Mt::ROWDIM, Mt::COLDIM, Mt>, Mt>::value)

// Ma and Mb are types
#define MATRIX_MUL_DIM_AGREE(Ma, Mb) \
            /* inner dimension match */ \
           (detail::_ImplMtxAdaptor<Ma, detail::ORIENT_VEC_ROW>::COLDIM == detail::_ImplMtxAdaptor<Mb, detail::ORIENT_VEC_COL>::ROWDIM or \
            /* ...or dynamic inner dimension demands runtime check: */ \
            detail::_ImplMtxAdaptor<Ma, detail::ORIENT_VEC_ROW>::COLDIM == DYNAMIC_DIM or \
            detail::_ImplMtxAdaptor<Mb, detail::ORIENT_VEC_COL>::ROWDIM == DYNAMIC_DIM)

// Ma and Mb are types
#define MATRIX_DIM_AGREE(Ma, Mb) \
    ((Ma::ROWDIM == Mb::ROWDIM or Ma::ROWDIM * Mb::ROWDIM == 0) and \
     (Ma::COLDIM == Mb::COLDIM or Ma::COLDIM * Mb::COLDIM == 0))

// Ma and Mb are types
#define LINALG_DIM_AGREE(Ma, Mb) \
       /* row dimension match */ \
    (((detail::_ImplMtxAdaptor<Ma, detail::_ImplVecOrient<Ma,Mb>::orient>::ROWDIM == detail::_ImplMtxAdaptor<Mb, detail::_ImplVecOrient<Ma,Mb>::orient>::ROWDIM) or \
      (detail::_ImplMtxAdaptor<Ma, detail::_ImplVecOrient<Ma,Mb>::orient>::ROWDIM  * detail::_ImplMtxAdaptor<Mb, detail::_ImplVecOrient<Ma,Mb>::orient>::ROWDIM == 0)) and \
       /* col dimension match */ \
     ((detail::_ImplMtxAdaptor<Ma, detail::_ImplVecOrient<Ma,Mb>::orient>::COLDIM == detail::_ImplMtxAdaptor<Mb, detail::_ImplVecOrient<Ma,Mb>::orient>::COLDIM) or \
      (detail::_ImplMtxAdaptor<Ma, detail::_ImplVecOrient<Ma,Mb>::orient>::COLDIM  * detail::_ImplMtxAdaptor<Mb, detail::_ImplVecOrient<Ma,Mb>::orient>::COLDIM == 0)))


#endif /* MATRIXMACRO_H_ */
