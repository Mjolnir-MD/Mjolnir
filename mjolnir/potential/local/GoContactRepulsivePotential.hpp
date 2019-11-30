#ifndef MJOLNIR_POTENTIAL_LOCAL_REPULSIVE_GO_CONTACT_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_REPULSIVE_GO_CONTACT_POTENTIAL_HPP

namespace mjolnir
{

template<typename T> class System;

// A repulsive portion of GoContactPotential.
// - V_repl(r) = e * (5*(r0/r)^12 - 6*(r0/r)^10 + 1) if r < r0
//             = 0                                   otherwise
//
template<typename realT>
class GoContactRepulsivePotential
{
  public:
    using real_type = realT;

  public:
    GoContactRepulsivePotential(const real_type e, const real_type v0) noexcept
        : epsilon_(e), v0_(v0)
    {}
    ~GoContactRepulsivePotential() = default;

    real_type potential(const real_type v) const noexcept
    {
        if(this->v0_ < v) {return real_type(0);}

        const auto vr  = this->v0_ / v;
        const auto vr3 = vr  * vr  * vr;
        const auto vr9 = vr3 * vr3 * vr3;
        return this->epsilon_ * (5 * vr9 * vr3 - 6 * vr9 * vr + real_type(1));
    }

    real_type derivative(const real_type v) const noexcept
    {
        if(this->v0_ < v) {return real_type(0);}

        const auto vinv = real_type(1) / v;
        const auto vr   = this->v0_ * vinv;
        const auto vr3  = vr  * vr  * vr;
        const auto vr9  = vr3 * vr3 * vr3;
        return this->epsilon_ * 60 * (vr9 * vr - vr9 * vr3) * vinv;
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "GoContactRepulsive";}

    real_type k()  const noexcept {return epsilon_;}
    real_type v0() const noexcept {return v0_;}

    real_type cutoff() const noexcept {return this->v0_;}

  private:

    real_type epsilon_;
    real_type v0_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GoContactRepulsivePotential<double>;
extern template class GoContactRepulsivePotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_POTENTIAL_LOCAL_GO_CONTACT_REPULTIVE_POTENTIAL */
