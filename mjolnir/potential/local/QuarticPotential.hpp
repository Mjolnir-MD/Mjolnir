#ifndef MJOLNIR_QUARTIC_POTENTIAL
#define MJOLNIR_QUARTIC_POTENTIAL

namespace mjolnir
{
template<typename T> class System;

/*! @brief quartic potential for 3SPN            *
 *  V(r) =     k1 * (r-r0)^2 +     k2 * (r-r0)^4 *
 * dV/dr = 2 * k1 * (r-r0)   + 4 * k2 * (r-r0)^3 */
template<typename traitsT>
class QuarticPotential
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    QuarticPotential(const real_type k1, const real_type k2,
                     const real_type native_val)
        : k1_(k1), k2_(k2), native_val_(native_val)
    {}
    ~QuarticPotential() = default;

    real_type potential(const real_type val) const
    {
        const real_type dval = val - this->native_val_;
        const real_type dv2 = dval * dval;
        return this->k1_ * dv2 + this->k2_ * dv2 * dv2;
    }

    real_type derivative(const real_type val) const
    {
        const real_type dr = val - this->native_val_;
        return  2 * this->k1_ * dr + 4 * this->k2_ * dr * dr * dr;
    }

    void update(const system_type&, const real_type) const noexcept {return;}

    const char* name() const noexcept {return "Quartic";}

  private:

    real_type k1_;
    real_type k2_;
    real_type native_val_; //!< most stable length
};

}
#endif /* MJOLNIR_QUARTIC_POTENTIAL */
