#ifndef MJOLNIR_HARMONIC_POTENTIAL
#define MJOLNIR_HARMONIC_POTENTIAL

namespace mjolnir
{

template<typename T> class System;

/*! @brief harmonic potential *
 *  V(r) = K * (r - r0)^2     *
 * dV/dr = 2 * K * (r - r0)   */
template<typename traitsT>
class HarmonicPotential
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

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

    void update(const system_type&, const real_type) const noexcept {return;}

    const char* name() const noexcept {return "Harmonic";}

  private:

    real_type k_;
    real_type r0_;
};

}
#endif /* MJOLNIR_HARMONIC_POTENTIAL */
