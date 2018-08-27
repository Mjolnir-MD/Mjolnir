#ifndef MJOLNIR_MATH_MATH_HPP
#define MJOLNIR_MATH_MATH_HPP
#include <cstdint>
#include <cmath>

namespace mjolnir
{

namespace compiletime
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

// after C++14, std::min and std::max become constexpr
template<typename T>
constexpr inline T min(const T x, const T y) noexcept
{
    return (x < y) ? x : y;
}

template<typename T>
constexpr inline T max(const T x, const T y) noexcept
{
    return (x < y) ? y : x;
}

// std::clamp is after C++17
template<typename T>
constexpr inline T clamp(const T x, const T low, const T high) noexcept
{
    return min(max(x, low), high);
}

} // compiletime

template<typename T>
inline T clamp(const T x, const T low, const T high) noexcept
{
    return std::min(std::max(x, low), high);
}

} // mjolnir
#endif // MJOLNIR_MATH_MATH_HPP
