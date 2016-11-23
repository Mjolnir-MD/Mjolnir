#ifndef MJOLNIR_POTENTIAL_BASE
#define MJOLNIR_POTENTIAL_BASE

namespace mjolnir
{

template<typename traitsT, std::size_t N>
struct PotentialBase
{
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordT;

    virtual real_type potential(std::array<coordT const&, N> pos) const = 0;
    virtual real_type derivative(std::array<coordT const&, N> pos) const = 0;
};

} // mjolnir

#endif /* MJOLNIR_POTENTIAL_BASE */
