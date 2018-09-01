#ifndef MJOLNIR_MATH_MATH_HPP
#define MJOLNIR_MATH_MATH_HPP
#include <mjolnir/math/compiletime.hpp>
#include <mjolnir/math/approx.hpp>
#include <mjolnir/math/Matrix.hpp>
#include <mjolnir/math/Vector.hpp>
#include <cstdint>
#include <cmath>

namespace mjolnir
{

// std::clamp is defined after C++17
// after C++14, std::min and std::max become constexpr
template<typename T>
inline T clamp(const T x, const T low, const T high) noexcept
{
    return std::min(std::max(x, low), high);
}

} // mjolnir
#endif // MJOLNIR_MATH_MATH_HPP
