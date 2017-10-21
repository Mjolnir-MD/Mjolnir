#ifndef MJOLNIR_MATH_VECTOR_H
#define MJOLNIR_MATH_VECTOR_H
#include <mjolnir/math/Matrix.hpp>
#include <mjolnir/math/quaternion.hpp>
#include <mjolnir/math/rsqrt.hpp>
#include <mjolnir/util/scalar_type_of.hpp>
#include <cmath>

namespace mjolnir
{

template<typename realT, std::size_t N>
using Vector = Matrix<realT, N, 1>;

// functions for vector 3d

template<typename realT>
inline realT
dot_product(const Vector<realT, 3>& lhs, const Vector<realT, 3>& rhs) noexcept
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

template<typename realT>
inline Vector<realT, 3>
cross_product(const Vector<realT, 3>& lhs, const Vector<realT, 3>& rhs) noexcept
{
    return Vector<realT, 3>(lhs[1] * rhs[2] - lhs[2] * rhs[1],
                            lhs[2] * rhs[0] - lhs[0] * rhs[2],
                            lhs[0] * rhs[1] - lhs[1] * rhs[0]);
}

template<typename realT>
inline realT
scalar_triple_product(const Vector<realT, 3>& lhs, const Vector<realT, 3>& mid,
                      const Vector<realT, 3>& rhs) noexcept
{
    return (lhs[1] * mid[2] - lhs[2] * mid[1]) * rhs[0] + 
           (lhs[2] * mid[0] - lhs[0] * mid[2]) * rhs[1] + 
           (lhs[0] * mid[1] - lhs[1] * mid[0]) * rhs[2];
}

template<typename realT>
inline realT length_sq(const Vector<realT, 3>& lhs) noexcept
{
    return lhs[0] * lhs[0] + lhs[1] * lhs[1] + lhs[2] * lhs[2];
}

template<typename realT>
inline realT length(const Vector<realT, 3>& lhs) noexcept
{
    return std::sqrt(length_sq(lhs));
}

template<typename realT>
Vector<realT>
rotate(const realT angle, const Vector<realT>& axis,
       const Vector<realT>& target) noexcept
{
    const realT half_angle     = angle / 2;
    const realT sin_normalized = std::sin(half_angle) * rsqrt(length_sq(axis));
    const quaternion<realT> Q{std::cos(half_angle),
                              axis[0] * sin_normalize,
                              axis[1] * sin_normalize,
                              axis[2] * sin_normalize};
    const quaternion<realT> P(0e0, target[0], target[1], target[2]);
    const quaternion<realT> S(Q * P * conj(Q));

    return Vector<realT>{S[1], S[2], S[3]};
}

}

#endif /* MJOLNIR_VECTOR */
