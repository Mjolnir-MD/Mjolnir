#ifndef MJOLNIR_POTENTIAL_LOCAL_ATTRACTIVE_GO_CONTACT_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_ATTRACTIVE_GO_CONTACT_POTENTIAL_HPP

namespace mjolnir
{

template<typename T> class System;

// An attractive portion of GoContactPotential.
// - V_attr(r) = -e                                            if r < r0_attr
//             = e * (5*(r0_attr/r)^12 - 6*(r0_attr/r)^10)     otherwise
//
template<typename realT>
class GoContactAttractivePotential
{
  public:
    using real_type = realT;

    static constexpr real_type cutoff_ratio()  {return 2.5;}

  public:
    GoContactAttractivePotential(const real_type e, const real_type v0) noexcept
        : epsilon_(e), v0_(v0)
    {}
    ~GoContactAttractivePotential() = default;

    real_type potential(const real_type v) const noexcept
    {
        if(v < this->v0_) {return -this->epsilon_;}
        if(this->cutoff() < v) {return real_type(0);}

        const auto vr  = this->v0_ / v;
        const auto vr3 = vr  * vr  * vr;
        const auto vr9 = vr3 * vr3 * vr3;
        return this->epsilon_ * (5 * vr9 * vr3 - 6 * vr9 * vr);
    }

    real_type derivative(const real_type v) const noexcept
    {
        if(v < this->v0_ || this->cutoff() < v) {return real_type(0);}

        const auto vinv = real_type(1) / v;
        const auto vr   = this->v0_ * vinv;
        const auto vr3  = vr  * vr  * vr;
        const auto vr9  = vr3 * vr3 * vr3;
        return this->epsilon_ * 60 * (vr9 * vr - vr9 * vr3) * vinv;
    }

    template<typename T>
    void initialize(const System<T>&) const noexcept {return;}

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "GoContactAttractive";}

    real_type k()  const noexcept {return epsilon_;}
    real_type v0() const noexcept {return v0_;}

    real_type cutoff()  const noexcept {return this->v0_ * cutoff_ratio();}

  private:

    real_type epsilon_;
    real_type v0_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GoContactAttractivePotential<double>;
extern template class GoContactAttractivePotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_POTENTIAL_LOCAL_ATTRACTIVE_GO_CONTACT_POTENTIAL */
