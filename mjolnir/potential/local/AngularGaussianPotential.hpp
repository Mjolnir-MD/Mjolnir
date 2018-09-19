#ifndef MJOLNIR_ANGULAR_GAUSSIAN_POTENTIAL
#define MJOLNIR_ANGULAR_GAUSSIAN_POTENTIAL
#include <mjolnir/math/math.hpp>

namespace mjolnir
{
template<typename T> class System;

// gaussian potential for dihedral interaction.
// V(r)  = epsilon * exp(-(r-r0)^2 / 2W^2)
// dV/dr = epsilon * (-(r-r0) / W^2) * exp(-(r-r0)^2 / 2W^2)
template<typename realT>
class AngularGaussianPotential
{
  public:
    using real_type = realT;

  public:
    AngularGaussianPotential(const real_type e, const real_type w,
                      const real_type v0) noexcept
        : epsilon_(e), inv_w2_(-1./(2*w*w)), v0_(v0)
    {}
    ~AngularGaussianPotential() = default;

    real_type potential(const real_type val) const noexcept
    {
        constexpr real_type pi     = math::constants<real_type>::pi;
        constexpr real_type two_pi = math::constants<real_type>::two_pi;

        real_type dval = val - this->v0_;
        if(dval < -pi) {dval += two_pi;} else if(pi < dval){dval -= two_pi;}
        return epsilon_ * std::exp(inv_w2_ * dval * dval);
    }

    real_type derivative(const real_type val) const noexcept
    {
        constexpr real_type pi     = math::constants<real_type>::pi;
        constexpr real_type two_pi = math::constants<real_type>::two_pi;

        real_type dval = val - this->v0_;
        if(dval < -pi) {dval += two_pi;} else if(pi < dval){dval -= two_pi;}
        return 2*inv_w2_ * dval * epsilon_ * std::exp(inv_w2_ * dval * dval);
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "AngularGaussian";}

  private:

    real_type epsilon_, inv_w2_;
    real_type v0_;
};

} // mjolnir
#endif /* MJOLNIR_GAUSSIAN_POTENTIAL */
