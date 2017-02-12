#ifndef MJOLNIR_VECTOR
#define MJOLNIR_VECTOR
#include "Matrix.hpp"
#include "quaternion.hpp"
#include "fast_inv_sqrt.hpp"
#include <mjolnir/util/scalar_type_of.hpp>
#include <cmath>

namespace mjolnir
{

template<typename realT, std::size_t N>
using Vector = Matrix<realT, N, 1>;

// for vector 3d
template<typename coordT>
inline typename scalar_type_of<coordT>::type
dot_product(const coordT& lhs, const coordT& rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

template<typename coordT>
inline coordT
cross_product(const coordT& lhs, const coordT& rhs)
{
    return coordT(lhs[1] * rhs[2] - lhs[2] * rhs[1],
                  lhs[2] * rhs[0] - lhs[0] * rhs[2],
                  lhs[0] * rhs[1] - lhs[1] * rhs[0]);
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
scalar_triple_product(const coordT& lhs, const coordT& mid, const coordT& rhs)
{
    return (lhs[1] * mid[2] - lhs[2] * mid[1]) * rhs[0] + 
           (lhs[2] * mid[0] - lhs[0] * mid[2]) * rhs[1] + 
           (lhs[0] * mid[1] - lhs[1] * mid[0]) * rhs[2];
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
length_sq(const coordT& lhs)
{
    return lhs[0] * lhs[0] + lhs[1] * lhs[1] + lhs[2] * lhs[2];
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
length(const coordT& lhs)
{
    return std::sqrt(length_sq(lhs));
}

template<typename coordT>
coordT
rotate(const typename scalar_type_of<coordT>::type angle,
       const coordT& axis, const coordT& target)
{
    typedef typename scalar_type_of<coordT>::type real_type;

    const real_type half_angle = angle * 0.5;
    const real_type sin_normalize =
        std::sin(half_angle) * fast_inv_sqrt(length_sq(axis));

    const quaternion<real_type> Q(std::cos(half_angle),
            axis[0] * sin_normalize, axis[1] * sin_normalize,
            axis[2] * sin_normalize);
    const quaternion<real_type> P(0e0, target[0], target[1], target[2]);
    const quaternion<real_type> S(Q * P * conj(Q));

    return coordT(S[1], S[2], S[3]);
}

}

#endif /* MJOLNIR_VECTOR */
