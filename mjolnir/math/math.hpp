#ifndef MJOLNIR_MATH_MATH_HPP
#define MJOLNIR_MATH_MATH_HPP
#include <cstdint>

// simple functions to calculate some parameters at a compile time

namespace mjolnir
{

template<typename realT>
constexpr inline realT pow(const realT x, const std::uint64_t N) noexcept
{
    return (N == 0) ? realT(1.0) : x * pow(x, N-1);
}

template<typename realT>
constexpr inline realT abs(const realT x) noexcept
{
    return (x > 0) ? x : -x;
}

} // mjolnir
#endif // MJOLNIR_MATH_MATH_HPP
