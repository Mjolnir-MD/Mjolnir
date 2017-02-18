#ifndef MJOLNIR_GLOBAL_POTENTIAL_BASE
#define MJOLNIR_GLOBAL_POTENTIAL_BASE
#include "ParticleContainer.hpp"
#include <string>

namespace mjolnir
{

template<typename traitsT>
class GlobalPotentialBase
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;

    virtual ~GlobalPotentialBase(){};

    virtual real_type potential(
        const std::size_t i, const std::size_t j, const real_type r) const = 0;
    virtual real_type derivative(
        const std::size_t i, const std::size_t j, const real_type r) const = 0;

    virtual real_type max_cutoff_length() const = 0;

    virtual void reset_parameter(const std::string&, const real_type) = 0;
};


} // mjolnir
#endif /* MJOLNIR_GLOBAL_POTENTIAL_BASE */
