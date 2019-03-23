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

template<typename charT, typename traitsT,
         typename T, std::size_t R>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const Matrix<T, R, 1>& mat)
{
    os << '(';
    for(std::size_t i=0; i<R; ++i)
    {
        if(i!=0) {os << ',';}
        os << mat(i, 0);
    }
    os << ')';
    return os;
}

// use mjolnir::math::X() to access elements of vector.

template<typename realT>
inline realT  X(Vector<realT, 3> const& v) noexcept {return v[0];}
template<typename realT>
inline realT& X(Vector<realT, 3>& v)       noexcept {return v[0];}

template<typename realT>
inline realT  Y(Vector<realT, 3> const& v) noexcept {return v[1];}
template<typename realT>
inline realT& Y(Vector<realT, 3>& v)       noexcept {return v[1];}

template<typename realT>
inline realT  Z(Vector<realT, 3> const& v) noexcept {return v[2];}
template<typename realT>
inline realT& Z(Vector<realT, 3>& v)       noexcept {return v[2];}

// functions for vector 3d

template<typename realT>
inline realT
dot_product(const Vector<realT, 3>& lhs, const Vector<realT, 3>& rhs) noexcept
{
    return X(lhs) * X(rhs) + Y(lhs) * Y(rhs) + Z(lhs) * Z(rhs);
}

template<typename realT>
inline Vector<realT, 3>
cross_product(const Vector<realT, 3>& lhs, const Vector<realT, 3>& rhs) noexcept
{
    return Vector<realT, 3>(Y(lhs) * Z(rhs) - Z(lhs) * Y(rhs),
                            Z(lhs) * X(rhs) - X(lhs) * Z(rhs),
                            X(lhs) * Y(rhs) - Y(lhs) * X(rhs));
}

template<typename realT>
inline realT length_sq(const Vector<realT, 3>& lhs) noexcept
{
    return X(lhs) * X(lhs) + Y(lhs) * Y(lhs) + Z(lhs) * Z(lhs);
}

template<typename realT>
inline realT length(const Vector<realT, 3>& lhs) noexcept
{
    return std::sqrt(length_sq(lhs));
}

} // math
} // mjolnir
#endif /* MJOLNIR_VECTOR */
