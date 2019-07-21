#ifndef MJOLNIR_POTENTIAL_GLOBAL_3SPN2_COMMON_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_3SPN2_COMMON_HPP
#include <cstdint>

namespace mjolnir
{

// X is for default parameter
enum class base_kind       : std::uint8_t {A,  T,  C,  G,  X};
enum class base_pair_kind  : std::uint8_t {AT, TA, CG, GC};
enum class base_stack_kind : std::uint8_t
{
    AA =  0, AT =  1, AC =  2, AG =  3,
    TA =  4, TT =  5, TC =  6, TG =  7,
    CA =  8, CT =  9, CC = 10, CG = 11,
    GA = 12, GT = 13, GC = 14, GG = 15,
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
    AA5 =  0, AT5 =  1, AC5 =  2, AG5 =  3,
    TA5 =  4, TT5 =  5, TC5 =  6, TG5 =  7,
    CA5 =  8, CT5 =  9, CC5 = 10, CG5 = 11,
    GA5 = 12, GT5 = 13, GC5 = 14, GG5 = 15,
    // sense 3 -- antisense 3 (Bj - Bi3)
    AA3 = 16, AT3 = 17, AC3 = 18, AG3 = 19,
    TA3 = 20, TT3 = 21, TC3 = 22, TG3 = 23,
    CA3 = 24, CT3 = 25, CC3 = 26, CG3 = 27,
    GA3 = 28, GT3 = 29, GC3 = 30, GG3 = 31
};

} // mjolnir
#endif// MJOLNIR_POTENTIAL_GLOBAL_3SPN2_COMMON_HPP
