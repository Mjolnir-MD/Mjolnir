#ifndef MJOLNIR_FORCEFIELD_3SPN2_COMMON_HPP
#define MJOLNIR_FORCEFIELD_3SPN2_COMMON_HPP
#include <iosfwd>
#include <limits>
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
enum class base_kind       : std::uint8_t {A=0, T=1, G=2, C=3, X};
enum class base_pair_kind  : std::uint8_t {AT, TA, GC, CG, INVALID};
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
    CA3 = 28, CT3 = 29, CG3 = 30, CC3 = 31,
    INVALID = 255
};

struct NucleotideInfo
{
    static constexpr std::size_t nil() noexcept
    {
        return std::numeric_limits<std::size_t>::max();
    }

    NucleotideInfo() noexcept
        : strand(nil()), P(nil()), S(nil()), B(nil()),
          base(base_kind::X)
    {}
    ~NucleotideInfo() noexcept = default;
    NucleotideInfo(NucleotideInfo const&) = default;
    NucleotideInfo(NucleotideInfo &&)     = default;
    NucleotideInfo& operator=(NucleotideInfo const&) = default;
    NucleotideInfo& operator=(NucleotideInfo &&)     = default;

    std::size_t strand;
    std::size_t P, S, B;
    base_kind   base;
};

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const bead_kind bk)
{
    switch(bk)
    {
        case bead_kind::Phosphate: {os << "Phosphate"; return os;}
        case bead_kind::Sugar:     {os << "Sugar"    ; return os;}
        case bead_kind::BaseA:     {os << "BaseA"    ; return os;}
        case bead_kind::BaseT:     {os << "BaseT"    ; return os;}
        case bead_kind::BaseG:     {os << "BaseG"    ; return os;}
        case bead_kind::BaseC:     {os << "BaseC"    ; return os;}
        case bead_kind::Unknown:   {os << "Unknown"  ; return os;}
        default:        {os << "invalid"  ; return os;}
    }
}
template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const base_kind bk)
{
    switch(bk)
    {
        case base_kind::A: {os << "A"; return os;}
        case base_kind::T: {os << "T"; return os;}
        case base_kind::G: {os << "G"; return os;}
        case base_kind::C: {os << "C"; return os;}
        case base_kind::X: {os << "X"; return os;}
        default:{os << "invalid"; return os;}
    }
}
template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const base_pair_kind bk)
{
    switch(bk)
    {
        case base_pair_kind::AT: {os << "AT"; return os;}
        case base_pair_kind::TA: {os << "TA"; return os;}
        case base_pair_kind::GC: {os << "GC"; return os;}
        case base_pair_kind::CG: {os << "CG"; return os;}
        default: {os << "invalid"; return os;}
    }
}
template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const base_stack_kind bs)
{
    switch(bs)
    {
        case base_stack_kind::AA: {os << "AA"; return os;}
        case base_stack_kind::AT: {os << "AT"; return os;}
        case base_stack_kind::AG: {os << "AG"; return os;}
        case base_stack_kind::AC: {os << "AC"; return os;}

        case base_stack_kind::TA: {os << "TA"; return os;}
        case base_stack_kind::TT: {os << "TT"; return os;}
        case base_stack_kind::TG: {os << "TG"; return os;}
        case base_stack_kind::TC: {os << "TC"; return os;}

        case base_stack_kind::GA: {os << "GA"; return os;}
        case base_stack_kind::GT: {os << "GT"; return os;}
        case base_stack_kind::GG: {os << "GG"; return os;}
        case base_stack_kind::GC: {os << "GC"; return os;}

        case base_stack_kind::CA: {os << "CA"; return os;}
        case base_stack_kind::CT: {os << "CT"; return os;}
        case base_stack_kind::CG: {os << "CG"; return os;}
        case base_stack_kind::CC: {os << "CC"; return os;}

        default: {os << "invalid"; return os;}
    }
}
template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const cross_stack_kind cs)
{
    switch(cs)
    {
        case cross_stack_kind::AA5: {os << "AA5"; return os;}
        case cross_stack_kind::AT5: {os << "AT5"; return os;}
        case cross_stack_kind::AG5: {os << "AG5"; return os;}
        case cross_stack_kind::AC5: {os << "AC5"; return os;}

        case cross_stack_kind::TA5: {os << "TA5"; return os;}
        case cross_stack_kind::TT5: {os << "TT5"; return os;}
        case cross_stack_kind::TG5: {os << "TG5"; return os;}
        case cross_stack_kind::TC5: {os << "TC5"; return os;}

        case cross_stack_kind::GA5: {os << "GA5"; return os;}
        case cross_stack_kind::GT5: {os << "GT5"; return os;}
        case cross_stack_kind::GG5: {os << "GG5"; return os;}
        case cross_stack_kind::GC5: {os << "GC5"; return os;}

        case cross_stack_kind::CA5: {os << "CA5"; return os;}
        case cross_stack_kind::CT5: {os << "CT5"; return os;}
        case cross_stack_kind::CG5: {os << "CG5"; return os;}
        case cross_stack_kind::CC5: {os << "CC5"; return os;}

        case cross_stack_kind::AA3: {os << "AA3"; return os;}
        case cross_stack_kind::AT3: {os << "AT3"; return os;}
        case cross_stack_kind::AG3: {os << "AG3"; return os;}
        case cross_stack_kind::AC3: {os << "AC3"; return os;}

        case cross_stack_kind::TA3: {os << "TA3"; return os;}
        case cross_stack_kind::TT3: {os << "TT3"; return os;}
        case cross_stack_kind::TG3: {os << "TG3"; return os;}
        case cross_stack_kind::TC3: {os << "TC3"; return os;}

        case cross_stack_kind::GA3: {os << "GA3"; return os;}
        case cross_stack_kind::GT3: {os << "GT3"; return os;}
        case cross_stack_kind::GG3: {os << "GG3"; return os;}
        case cross_stack_kind::GC3: {os << "GC3"; return os;}

        case cross_stack_kind::CA3: {os << "CA3"; return os;}
        case cross_stack_kind::CT3: {os << "CT3"; return os;}
        case cross_stack_kind::CG3: {os << "CG3"; return os;}
        case cross_stack_kind::CC3: {os << "CC3"; return os;}

        default: {os << "invalid"; return os;}
    }
}

template<typename charT, typename traitsT>
std::basic_ostream<charT, traitsT>&
operator<<(std::basic_ostream<charT, traitsT>& os, const NucleotideInfo& ni)
{
    os << "{strand = " << ni.strand
       << ", P = " << ni.P << ", S = " << ni.S << ", B = " << ni.B << "}";
    return os;
}

} // 3SPN2
} // mjolnir
#endif// MJOLNIR_POTENTIAL_GLOBAL_3SPN2_COMMON_HPP
