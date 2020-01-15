#ifndef MJOLNIR_UTIL_IS_FINITE_HPP
#define MJOLNIR_UTIL_IS_FINITE_HPP
#include <cstring>
#include <cstdint>

namespace mjolnir
{

// When we compile a simulator, we often turn optimization flags on.
// Some flags that allow relatively agressive optimizations, e.g. -Ofast and
// -ffast-math, NaN-cheking will be elimited. It is because those flags
// implicitly turns a flag `-ffinite-math-only`, at least on GCC, that assumes
// that all the floating point value will never be NaN. Under that condition,
// naive NaN-checking, such as calling std::isnan, will be eliminated as a
// dead code. But sometimes we want to check it anyway.
//     To avoid that optimization, those `mjolnir::is_finite` functions directly
// check bits of a value. By definition, all the exponent bits of NaN and Inf
// are 1. It checks those bits by bit operations. It can take longer time than
// __builtin_isnan that compares 2 floating point numbers. So don't use this
// in hotspots.

inline bool is_finite(const double x) noexcept
{
    std::uint64_t n;
    std::memcpy(reinterpret_cast<char*>(&n),
                reinterpret_cast<const char*>(&x), sizeof(double));

    // bin: 0111'1111'1111'0000'0000'...'0000
    // hex:    7    F    F    0    0 ...    0
    constexpr std::uint64_t mask = 0x7FF0000000000000;
    return ((n & mask) ^ mask) != 0u;
}

inline bool is_finite(const float x) noexcept
{
    std::uint32_t n;
    std::memcpy(reinterpret_cast<char*>(&n),
                reinterpret_cast<const char*>(&x), sizeof(float));

    // bin: 0111'1111'1000'0000'0000'...'0000
    // hex:    7    F    8    0    0 ...    0
    constexpr std::uint32_t mask = 0x7F800000;
    return ((n & mask) ^ mask) != 0u;
}
} // mjolnir
#endif // MJOLNIR_UTIL_IS_FINITE_HPP
