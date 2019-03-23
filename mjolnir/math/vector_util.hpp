#ifndef MJOLNIR_MATH_VECTOR_UTIL_H
#define MJOLNIR_MATH_VECTOR_UTIL_H

// these codes are for buffering difference between implementations.

namespace mjolnir
{
namespace math
{

template<typename CoordinateT>
struct real_type_of;

template<typename CoordinateT>
using real_type_of_t = typename real_type_of<CoordinateT>::type;

template<typename CoordinateT>
struct make_coordinate_impl;

template<typename CoordinateT>
inline CoordinateT make_coordinate(const real_type_of_t<CoordinateT> x,
                                   const real_type_of_t<CoordinateT> y,
                                   const real_type_of_t<CoordinateT> z) noexcept
{
    return make_coordinate_impl<CoordinateT>::invoke(x, y, z);
}

} // math
} // mjolnir
#endif//MJOLNIR_MATH_VECTOR_UTIL_H
