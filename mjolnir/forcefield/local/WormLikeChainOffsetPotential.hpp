#ifndef MJOLNIR_POTENTIAL_LOCAL_WORM_LIKE_CHAIN_OFFSET_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_LOCAL_WORM_LIKE_CHAIN_OFFSET_POTENTIAL_HPP
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

// Worm-Like chain type interaction, which is usually applyed to the interaction
// via polypeptide linker. The formula of this potential partialy referenced
// Guillaume et al., (2013).
// This interaction include offset for r.
// V(r)  = kT/p  * (lc/4 * (1/(1 - (r - offset)/lc) - 1) - (r - offset)/4 + (r - offset)^2/2lc)
// dV/dr = kT/p * (1/4 * (1/(1 - (r - offset)/lc)^2 - 1) - (r - offset)/lc)
template<typename realT>
class WormLikeChainOffsetPotential
{
  public:
    using real_type   = realT;

  public:
    WormLikeChainOffsetPotential(const real_type p, const real_type lc, const real_type offset) noexcept
        : p_(p), lc_(lc), inv_lc_(1.0 / lc), offset_(offset), temperature_(300.0)
    {
        // XXX should be updated before use because T is default value!
        this->calc_parameters();
    }
    ~WormLikeChainOffsetPotential() = default;

    real_type potential(const real_type dist) const noexcept
    {
        if(dist <= offset_) {return kBT_4p_ * lc_;}
        const real_type l = dist - offset_;
        return kBT_4p_ * (lc_* (lc_ / (lc_ - l) - 1.0) - l + 2.0 * l * l * inv_lc_);
    }

    real_type derivative(const real_type dist) const noexcept
    {
        if(dist <= offset_) {return 0.0;}
        const real_type l      = dist - offset_;
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
            MJOLNIR_LOG_ERROR("WormLikeChainOffset requires `tempereture` attribute");
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

    static const char* name() noexcept {return "WormLikeChainOffset";}

    real_type lc() const noexcept {return lc_;}
    real_type p() const noexcept {return p_;}
    real_type offset() const noexcept {return offset_;}

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
    real_type offset_;
    real_type temperature_; // [K]
    real_type kBT_4p_; // kBT / 4p
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class WormLikeChainOffsetPotential<double>;
extern template class WormLikeChainOffsetPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
# endif /* MJOLNIR_POTENTIAL_LOCAL_WORM_LIKE_CHAIN_OFFSET_POTENTIAL_HPP */
