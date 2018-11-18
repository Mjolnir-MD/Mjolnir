#ifndef MJOLNIR_HARMONIC_POTENTIAL
#define MJOLNIR_HARMONIC_POTENTIAL

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
    HarmonicPotential(const real_type k, const real_type r0)
        : k_(k), r0_(r0)
    {}
    ~HarmonicPotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        const real_type dr = r - r0_;
        return this->k_ * dr * dr;
    }

    real_type derivative(const real_type r) const noexcept
    {
        return  2 * this->k_ * (r - r0_);
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Harmonic";}

  private:

    real_type k_;
    real_type r0_;
};

}
#endif /* MJOLNIR_HARMONIC_POTENTIAL */
