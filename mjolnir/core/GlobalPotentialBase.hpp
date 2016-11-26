#ifndef MJOLNIR_GLOBAL_POTENTIAL_BASE
#define MJOLNIR_GLOBAL_POTENTIAL_BASE

namespace mjolnir
{

template<typename traitsT, typename valueT>
class GlobalPotentialBase
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;

    virtual ~GlobalPotentialBase(){};

    virtual real_type
    potential(const valueT val1, const valueT va2, const real_type dist) = 0;
    virtual real_type
    derivative(const valueT val1, const valueT va2, const real_type dist) = 0;
};


} // mjolnir
#endif /* MJOLNIR_GLOBAL_POTENTIAL_BASE */
