#ifndef MJOLNIR_POTENTIAL_LOCAL_ANGULAR_GAUSSIAN_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_ANGULAR_GAUSSIAN_POTENTIAL_HPP
#include <mjolnir/math/math.hpp>
#include <limits>

namespace mjolnir
{
template<typename T> class System;

// gaussian potential for dihedral interaction.
// V(r)  = epsilon * exp(-(r-r0)^2 / 2W^2)
// dV/dr = epsilon * (-(r-r0) / W^2) * exp(-(r-r0)^2 / 2W^2)
template<typename realT>
class PeriodicGaussianPotential
{
  public:
    using real_type = realT;

  public:
    PeriodicGaussianPotential(const real_type k, const real_type sigma,
                      const real_type v0) noexcept
        : k_(k), sigma_(sigma), inv_sigma2_(-1./(2*sigma*sigma)), v0_(v0)
    {}
    ~PeriodicGaussianPotential() = default;

    real_type potential(const real_type val) const noexcept
    {
        constexpr real_type pi     = math::constants<real_type>::pi();
        constexpr real_type two_pi = math::constants<real_type>::two_pi();

        real_type dval = val - this->v0_;
        if(dval < -pi) {dval += two_pi;} else if(pi < dval){dval -= two_pi;}
        return k_ * std::exp(inv_sigma2_ * dval * dval);
    }

    real_type derivative(const real_type val) const noexcept
    {
        constexpr real_type pi     = math::constants<real_type>::pi();
        constexpr real_type two_pi = math::constants<real_type>::two_pi();

        real_type dval = val - this->v0_;
        if(dval < -pi) {dval += two_pi;} else if(pi < dval){dval -= two_pi;}
        return 2*inv_sigma2_ * dval * k_ * std::exp(inv_sigma2_ * dval * dval);
    }

    template<typename T>
    void initialize(const System<T>&) const noexcept {return;}

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "PeriodicGaussian";}

    real_type k()     const noexcept {return k_;}
    real_type sigma() const noexcept {return sigma_;}
    real_type v0()    const noexcept {return v0_;}

    real_type cutoff() const noexcept // no cutoff exists.
    {return std::numeric_limits<real_type>::infinity();}

  private:

    real_type k_, sigma_;
    real_type inv_sigma2_;
    real_type v0_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class PeriodicGaussianPotential<double>;
extern template class PeriodicGaussianPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_GAUSSIAN_POTENTIAL */
