#ifndef MJOLNIR_MATH_VECTOR_HPP
#define MJOLNIR_MATH_VECTOR_HPP
#include <mjolnir/math/Matrix.hpp>
#include <mjolnir/math/functions.hpp>
#include <mjolnir/math/vector_util.hpp>
#include <cmath>

namespace mjolnir
{
namespace math
{

template<typename realT, std::size_t N>
using Vector = Matrix<realT, N, 1>;

template<typename realT, std::size_t N>
struct real_type_of<Vector<realT, N>>
{
    using type = realT;
};

template<typename realT, std::size_t N>
struct make_coordinate_impl<Vector<realT, N>>
{
    static Vector<realT, N>
    invoke(const realT x, const realT y, const realT z) noexcept
    {
        return Vector<realT, N>(x, y, z);
    }
};

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

template<typename charT, typename traitsT, typename T>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const Matrix<T, 3, 1>& vec)
{
    os << '(' << X(vec) << ',' << Y(vec) << ',' << Z(vec) << ')';
    return os;
}

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

template<typename realT>
inline realT rlength(const Vector<realT, 3>& lhs) noexcept
{
    return ::mjolnir::math::rsqrt(length_sq(lhs));
}

template<typename realT>
inline Matrix<realT, 3, 3> tensor_product(
        const Vector<realT, 3>& lhs, const Vector<realT, 3>& rhs) noexcept
{
    return Matrix<realT, 3, 3>(
            X(lhs) * X(rhs), X(lhs) * Y(rhs), X(lhs) * Z(rhs),
            Y(lhs) * X(rhs), Y(lhs) * Y(rhs), Y(lhs) * Z(rhs),
            Z(lhs) * X(rhs), Z(lhs) * Y(rhs), Z(lhs) * Z(rhs));
}

} // math
} // mjolnir
#endif /* MJOLNIR_VECTOR */
