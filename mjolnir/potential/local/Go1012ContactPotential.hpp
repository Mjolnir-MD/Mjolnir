#ifndef MJOLNIR_GO_10_12_CONTACT_POTENTIAL
#define MJOLNIR_GO_10_12_CONTACT_POTENTIAL

namespace mjolnir
{

template<typename T> class System;

/*! @brief Go contact 10-12 potential *
 *  V(r) = epsilon * (5 * (r0/r)^12 - 6 * (r0/r)^10)  *
 * dV/dr = 60 * epsilon * ((r0/r)^10 - (r0/r)^12) / r */
template<typename realT>
class Go1012ContactPotential
{
  public:
    using real_type = realT;

    static constexpr real_type cutoff_ratio  = 2.5;
    static constexpr real_type rcutoff_ratio = 1. / cutoff_ratio;

  public:
    Go1012ContactPotential(const real_type e, const real_type r0) noexcept
        : epsilon_(e), r0_(r0)
    {}
    ~Go1012ContactPotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        const real_type rd   = this->r0_ / r;
        if(rd < rcutoff_ratio) return 0.0;

        const real_type rd2  = rd  * rd;
        const real_type rd4  = rd2 * rd2;
        const real_type rd8  = rd4 * rd4;
        const real_type rd10 = rd8 * rd2;
        const real_type rd12 = rd8 * rd4;

        return this->epsilon_ * (5. * rd12 - 6. * rd10);
    }

    real_type derivative(const real_type r) const noexcept
    {
        const real_type invr = 1. / r;
        const real_type rd   = this->r0_ * invr;
        if(rd < rcutoff_ratio) return 0.0;

        const real_type rd2  = rd  * rd;
        const real_type rd4  = rd2 * rd2;
        const real_type rd8  = rd4 * rd4;
        const real_type rd10 = rd8 * rd2;
        const real_type rd12 = rd8 * rd4;
        return this->epsilon_ * 60. * (rd10 - rd12) * invr;
    }

    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    static const char* name() noexcept {return "Go1012Contact";}

  private:

    real_type epsilon_;
    real_type r0_;
};

template<typename realT>
constexpr typename Go1012ContactPotential<realT>::real_type
Go1012ContactPotential<realT>::cutoff_ratio;

template<typename realT>
constexpr typename Go1012ContactPotential<realT>::real_type
Go1012ContactPotential<realT>::rcutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_GO_10_12_CONTACT_POTENTIAL */
