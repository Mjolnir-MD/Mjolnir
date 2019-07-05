#ifndef MJOLNIR_POTENTIAL_LOCAL_GAUSSIAN_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_GAUSSIAN_POTENTIAL_HPP
#include <mjolnir/math/constants.hpp>
#include <cmath>

namespace mjolnir
{
template<typename T> class System;

// well-known gaussian-shape potential.
// V(r)  = k * exp(-(r-r0)^2 / 2sigma^2)
// dV/dr = k * (-(r-r0) / sigma^2) * exp(-(r-r0)^2 / 2sigma^2)
template<typename realT>
class GaussianPotential
{
  public:
    using real_type = realT;

  public:
    GaussianPotential(const real_type k, const real_type sigma,
                      const real_type v0) noexcept
        : k_(k), sigma_(sigma), inv_sigma2_(-1./(2*sigma*sigma)), v0_(v0),
          cutoff_(v0_ + sigma_ *
                  std::sqrt(2 * std::log(k_ / math::tolerance<real_type>())))
    {}
    ~GaussianPotential() = default;

    real_type potential(const real_type val) const noexcept
    {
        const real_type dval = val - this->v0_;
        return k_ * std::exp(inv_sigma2_ * dval * dval);
    }

    real_type derivative(const real_type val) const noexcept
    {
        const real_type dval = val - this->v0_;
        return 2*inv_sigma2_ * dval * k_ * std::exp(inv_sigma2_ * dval * dval);
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Gaussian";}

    real_type k()      const noexcept {return k_;}
    real_type sigma()  const noexcept {return sigma_;}
    real_type v0()     const noexcept {return v0_;}

    real_type cutoff() const noexcept {return cutoff_;}

  private:

    real_type k_, sigma_;
    real_type inv_sigma2_;
    real_type v0_;
    real_type cutoff_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GaussianPotential<double>;
extern template class GaussianPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_GAUSSIAN_POTENTIAL */
