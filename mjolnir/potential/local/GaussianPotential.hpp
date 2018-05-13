#ifndef MJOLNIR_GAUSSIAN_POTENTIAL
#define MJOLNIR_GAUSSIAN_POTENTIAL
#include <cmath>

namespace mjolnir
{
template<typename T> class System;

/*! @brief gaussian potential for AICG2+                     *
 * V(r)  = epsilon * exp(-(r-r0)^2 / 2W^2)                   *
 * dV/dr = epsilon * (-(r-r0) / W^2) * exp(-(r-r0)^2 / 2W^2) */
template<typename traitsT>
class GaussianPotential
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    GaussianPotential(const real_type e, const real_type w,
                      const real_type eq_val) noexcept
        : epsilon_(e), inv_w2_(-1./(2*w*w)), eq_val_(eq_val)
    {}
    ~GaussianPotential() = default;

    real_type potential(const real_type val) const noexcept
    {
        const real_type dval = val - this->eq_val_;
        return epsilon_ * std::exp(inv_w2_ * dval * dval);
    }

    real_type derivative(const real_type val) const noexcept
    {
        const real_type dval = val - this->eq_val_;
        return 2*inv_w2_ * dval * epsilon_ * std::exp(inv_w2_ * dval * dval);
    }

    void update(const system_type&, const real_type) const noexcept {return;}

    const char* name() const noexcept {return "Gaussian";}

  private:

    real_type epsilon_, inv_w2_;
    real_type eq_val_;
};

} // mjolnir
#endif /* MJOLNIR_GAUSSIAN_POTENTIAL */
