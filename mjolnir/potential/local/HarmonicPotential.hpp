#ifndef MJOLNIR_POTENTIAL_LOCAL_HARMONIC_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_HARMONIC_POTENTIAL_HPP
#include <limits>

namespace mjolnir
{

template<typename T> class System;

// Well-known harmonic potential.
template<typename realT>
class HarmonicPotential
{
  public:
    using real_type = realT;

  public:
    HarmonicPotential(const real_type k, const real_type v0)
        : k_(k), v0_(v0)
    {}
    ~HarmonicPotential() = default;

    real_type potential(const real_type v) const noexcept
    {
        const real_type dv = v - v0_;
        return this->k_ * dv * dv;
    }

    real_type derivative(const real_type v) const noexcept
    {
        return  2 * this->k_ * (v - v0_);
    }

    template<typename T>
    void initialize(const System<T>&) const noexcept {return;}

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Harmonic";}

    real_type k()  const noexcept {return k_;}
    real_type v0() const noexcept {return v0_;}

    real_type cutoff() const noexcept // no cutoff exists.
    {return std::numeric_limits<real_type>::infinity();}

  private:

    real_type k_;
    real_type v0_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class HarmonicPotential<double>;
extern template class HarmonicPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

}
#endif /* MJOLNIR_HARMONIC_POTENTIAL */
