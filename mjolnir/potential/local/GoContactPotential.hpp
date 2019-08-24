#ifndef MJOLNIR_POTENTIAL_LOCAL_GO_CONTACT_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_GO_CONTACT_POTENTIAL_HPP

namespace mjolnir
{

template<typename T> class System;

// Go-like contact potential that is a kind of 10-12 Lennard-Jones interaction
// used in off-lattice Go-like protein model (Clement et al., 2000).
//  V(r) = epsilon * (5 * (r0/r)^12 - 6 * (r0/r)^10)
// dV/dr = 60 * epsilon * ((r0/r)^10 - (r0/r)^12) / r
template<typename realT>
class GoContactPotential
{
  public:
    using real_type = realT;
    static constexpr real_type cutoff_ratio() {return 2.5;}

  public:
    GoContactPotential(const real_type e, const real_type v0) noexcept
        : epsilon_(e), v0_(v0)
    {}
    ~GoContactPotential() = default;

    real_type potential(const real_type v) const noexcept
    {
        if(this->cutoff() <= v) {return real_type(0.0);}

        const real_type rd   = this->v0_ / v;
        const real_type rd3  = rd  * rd  * rd;
        const real_type rd9  = rd3 * rd3 * rd3;

        return this->epsilon_ * (5 * rd9 * rd3 - 6 * rd9 * rd);
    }

    real_type derivative(const real_type v) const noexcept
    {
        if(this->cutoff() <= v) {return real_type(0.0);}

        const real_type invr = real_type(1.0) / v;
        const real_type rd   = this->v0_ * invr;
        const real_type rd3  = rd  * rd  * rd;
        const real_type rd9  = rd3 * rd3 * rd3;

        return this->epsilon_ * 60 * (rd9 * rd - rd9 * rd3) * invr;
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "GoContact";}

    real_type k()      const noexcept {return epsilon_;}
    real_type v0()     const noexcept {return v0_;}

    real_type cutoff() const noexcept {return this->v0_ * cutoff_ratio();}

  private:

    real_type epsilon_;
    real_type v0_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GoContactPotential<double>;
extern template class GoContactPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_GO_10_12_CONTACT_POTENTIAL */
