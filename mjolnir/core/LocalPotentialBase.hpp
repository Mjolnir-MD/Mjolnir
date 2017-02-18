#ifndef MJOLNIR_LOCAL_POTENTIAL_BASE
#define MJOLNIR_LOCAL_POTENTIAL_BASE
#include <string>

namespace mjolnir
{

template<typename traitsT>
class LocalPotentialBase
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;

  public:
    LocalPotentialBase() = default;
    virtual ~LocalPotentialBase() = default;

    virtual real_type potential(const real_type value) const = 0;
    virtual real_type derivative(const real_type value) const = 0;

    virtual void reset_parameter(const std::string&, const real_type) = 0;
};

}// mjolnir
#endif /* MJOLNIR_LOCAL_POTENTIAL_BASE */
