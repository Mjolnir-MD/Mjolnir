#ifndef MJOLNIR_HARMONIC_POTENTIAL
#define MJOLNIR_HARMONIC_POTENTIAL
#include "LocalPotentialBase.hpp"

namespace mjolnir
{

/*! @brief harmonic potential & derivative                              *
 * designed for local force field. so it has its own parameters.        *
 * DO NOT USE THIS FOR GLOBAL FORCE FIELD! this may occupy your memory. *
 * V(phi)  = K * (value - native_value)^2                               *
 * dV/dphi = 2 * K * (value - native_value)                             */
template<typename traitsT>
class HarmonicPotential : public LocalPotentialBase<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    HarmonicPotential(const real_type k, const real_type native_val)
        : k_(k), native_val_(native_val)
    {}
    ~HarmonicPotential() override = default;

    real_type potential(const real_type val) const override
    {
        const real_type dval = val - this->native_val_;
        return this->k_ * dval * dval;
    }

    real_type derivative(const real_type val) const override
    {
        return  2 * this->k_ * (val - this->native_val_);
    }

  private:

    const real_type k_;          //!< spring coefficient
    const real_type native_val_; //!< most stable length
};

}
#endif /* MJOLNIR_HARMONIC_POTENTIAL */
