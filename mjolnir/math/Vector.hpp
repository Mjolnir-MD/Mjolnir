#ifndef MJOLNIR_MATH_VECTOR_H
#define MJOLNIR_MATH_VECTOR_H
#include <mjolnir/math/Matrix.hpp>
#include <mjolnir/math/quaternion.hpp>
#include <mjolnir/math/functions.hpp>
#include <cmath>

namespace mjolnir
{
namespace math
{

template<typename realT, std::size_t N>
using Vector = Matrix<realT, N, 1>;

// use mjolnir::math::X() to access elements of vector.

template<typename realT>
inline realT  X(const Vector<realT, 3>& v) noexcept {return v[0];}
template<typename realT>
inline realT& X(Vector<realT, 3>& v) noexcept {return v[0];}

template<typename realT>
inline realT  Y(const Vector<realT, 3>& v) noexcept {return v[1];}
template<typename realT>
inline realT& Y(Vector<realT, 3>& v) noexcept {return v[1];}

template<typename realT>
inline realT  Z(const Vector<realT, 3>& v) noexcept {return v[2];}
template<typename realT>
inline realT& Z(Vector<realT, 3>& v) noexcept {return v[2];}

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

} // math
} // mjolnir
#endif /* MJOLNIR_VECTOR */
