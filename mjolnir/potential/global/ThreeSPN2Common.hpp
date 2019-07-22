#ifndef MJOLNIR_POTENTIAL_GLOBAL_3SPN2_COMMON_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_3SPN2_COMMON_HPP
#include <cstdint>

namespace mjolnir
{
namespace parameter_3SPN2
{
enum class bead_kind: std::uint8_t
{
    Phosphate, Sugar, BaseA, BaseT, BaseG, BaseC, Unknown
};

// X is for default parameter
enum class base_kind       : std::uint8_t {A=0, T=1, G=2, C=3,  X};
enum class base_pair_kind  : std::uint8_t {AT, TA, GC, CG};
enum class base_stack_kind : std::uint8_t
{
    AA =  0, AT =  1, AG =  2, AC =  3,
    TA =  4, TT =  5, TG =  6, TC =  7,
    GA =  8, GT =  9, GG = 10, GC = 11,
    CA = 12, CT = 13, CG = 14, CC = 15
};
enum class cross_stack_kind: std::uint8_t
{
    //       Si   Bi   Bj   Sj
    //  5'    o -- o===o -- o     3'
    //  ^    /      \ /      \    |
    //  | P o        x        o P |
    //  |    \      / \      /    v
    //  3'    o -- o===o -- o     5'
    //           Bi3   Bj5
    //
    // sense 5 -- antisense 5 (Bi - Bj5)
    AA5 =  0, AT5 =  1, AG5 =  2, AC5 =  3,
    TA5 =  4, TT5 =  5, TG5 =  6, TC5 =  7,
    GA5 =  8, GT5 =  9, GG5 = 10, GC5 = 11,
    CA5 = 12, CT5 = 13, CG5 = 14, CC5 = 15,
    // sense 3 -- antisense 3 (Bj - Bi3)
    AA3 = 16, AT3 = 17, AG3 = 18, AC3 = 19,
    TA3 = 20, TT3 = 21, TG3 = 22, TC3 = 23,
    GA3 = 24, GT3 = 25, GG3 = 26, GC3 = 27,
    CA3 = 28, CT3 = 29, CG3 = 30, CC3 = 31
};

} // 3SPN2
} // mjolnir
#endif// MJOLNIR_POTENTIAL_GLOBAL_3SPN2_COMMON_HPP
