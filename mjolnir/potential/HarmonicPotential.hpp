#ifndef MJOLNIR_HARMONIC_POTENTIAL
#define MJOLNIR_HARMONIC_POTENTIAL
#include <mjolnir/core/LocalPotentialBase.hpp>

namespace mjolnir
{

/*! @brief harmonic potential *
 *  V(r) = K * (r - r0)^2     *
 * dV/dr = 2 * K * (r - r0)   */
template<typename traitsT>
class HarmonicPotential : public LocalPotentialBase<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    HarmonicPotential(const real_type k, const real_type r0)
        : k_(k), r0_(r0)
    {}
    ~HarmonicPotential() override = default;

    real_type potential(const real_type r) const override
    {
        const real_type dr = r - r0_;
        return this->k_ * dr * dr;
    }

    real_type derivative(const real_type r) const override
    {
        return  2 * this->k_ * (r - r0_);
    }

  private:

    const real_type k_;
    const real_type r0_;
};

}
#endif /* MJOLNIR_HARMONIC_POTENTIAL */
