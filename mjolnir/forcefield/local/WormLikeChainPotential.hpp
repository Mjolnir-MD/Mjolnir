#ifndef MJOLNIR_POTENTIAL_LOCAL_WORM_LIKE_CHAIN_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_WORM_LIKE_CHAIN_POTENTIAL_HPP
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

// Worm-Like chain type interaction, which is usually applyed to the interaction
// via polypeptide linker. The formula of this potential partialy referenced
// Guillaume et al., (2013).
// V(r)  = kT/p  * (lc/4 * 1/(1 - r/rc) - r/4 + r^2/2lc)
// dV/dr = kT/p * (1/4 * (1/(1 - r/lc)^2 - 1) - r/lc)
template<typename realT>
class WormLikeChainPotential
{
  public:
    using real_type   = realT;

  public:
    WormLikeChainPotential(const real_type p, const real_type lc) noexcept
        : p_(p), lc_(lc), inv_lc_(1.0 / lc), temperature_(300.0)
    {
        // XXX should be updated before use because T is default value!
        this->calc_parameters();
    }
    ~WormLikeChainPotential() = default;

    real_type potential(const real_type l) const noexcept
    {
        return kBT_4p_ * (1.0 / (lc_ - l) - lc_ - l + 2.0 * l * l * inv_lc_);
    }

    real_type derivative(const real_type l) const noexcept
    {
        const real_type l_lc   = l * inv_lc_;

        const real_type diff_one_l_lc = 1 - l_lc;
        return kBT_4p_ * (1.0 / (diff_one_l_lc * diff_one_l_lc) - 1.0 + 4.0 * l_lc);
    }

    template<typename T>
    void initialize(const System<T>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(!sys.has_attribute("temperature"))
        {
            MJOLNIR_LOG_ERROR("WormLikeChain requires `tempereture` attribute");
        }
        this->update(sys); // calc parameters
        return;
    }

    // for temperature chainges...
    template<typename T>
    void update(const System<T>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        assert(sys.has_attribute("temperature"));
        temperature_ = sys.attribute("temperature");

        MJOLNIR_LOG_INFO("temperature = ", temperature_);

        calc_parameters();
        return;
    }

    static const char* name() noexcept {return "WormLikeChain";}

    real_type lc() const noexcept {return lc_;}
    real_type p() const noexcept {return p_;}

    real_type cutoff() const noexcept // no cutoff exists.
    {return std::numeric_limits<real_type>::infinity();}

  private:

    void calc_parameters() noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        const real_type kB = physics::constants<real_type>::kB();

        MJOLNIR_LOG_INFO("kB = ", kB);

        kBT_4p_ = (kB * temperature_) / (4.0 * p_);
        return;
    }

  private:

    real_type p_, lc_;
    real_type inv_lc_; // 1 / lc
    real_type temperature_; // [K]
    real_type kBT_4p_; // kBT / 4p
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class WormLikeChainPotential<double>;
extern template class WormLikeChainPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
# endif /* MJOLNIR_POTENTIAL_LOCAL_WORM_LIKE_CHAIN_POTENTIAL_HPP */
