#ifndef MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#define MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/potential/global/ChainIgnoration.hpp>
#include <mjolnir/math/constants.hpp>
#include <vector>
#include <cmath>

namespace mjolnir
{

template<typename realT, typename ChainIgnoration>
class DebyeHuckelPotential
{
  public:
    using real_type             = realT;
    using container_type        = std::vector<real_type>;
    using chain_ignoration_type = ChainIgnoration;
    using topology_type         = Topology;
    using chain_id_type         = typename topology_type::chain_id_type;
    using connection_kind_type  = typename topology_type::connection_kind_type;

    // r_cutoff = cutoff_ratio * debye_length
    static constexpr real_type cutoff_ratio = 5.5;

  public:

    DebyeHuckelPotential(const container_type& charges,
        const std::map<connection_kind_type, std::size_t>& exclusions)
        : charges_(charges), temperature_(300.0), ion_conc_(0.1),
          ignore_within_(exclusions.begin(), exclusions.end())
    {
        // XXX should be updated before use because T and ion conc are default!
        this->calc_parameters();
    }
    DebyeHuckelPotential(container_type&& charges,
        const std::map<connection_kind_type, std::size_t>& exclusions)
        : charges_(std::move(charges)), temperature_(300.0), ion_conc_(0.1),
          ignore_within_(exclusions.begin(), exclusions.end())
    {
        // XXX should be updated before use because T and ion conc are default!
        this->calc_parameters();
    }
    ~DebyeHuckelPotential() = default;

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const noexcept
    {
        if(this->max_cutoff_length() <= r) {return 0.0;}
        return this->inv_4_pi_eps0_epsk_ *
               this->charges_[i] * this->charges_[j] *
               std::exp(-r * this->inv_debye_length_) / r;
    }

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        if(this->max_cutoff_length() <= r) {return 0.0;}
        return -(this->inv_4_pi_eps0_epsk_) *
               (debye_length_ + r) * this->inv_debye_length_ / (r * r) *
               this->charges_[i] * this->charges_[j] *
               std::exp(-r * this->inv_debye_length_);
    }

    real_type max_cutoff_length() const noexcept
    {
        return debye_length_ * cutoff_ratio;
    }

    // for temperature/ionic concentration changes...
    template<typename traitsT>
    void update(const System<traitsT>& sys) noexcept
    {
        this->temperature_ = sys.attribute("temperature");
        this->ion_conc_    = sys.attribute("ionic_strength");
        this->calc_parameters();
        return;
    }

    // e.g. {"bond", 3} means ignore particles connected within 3 "bond"s
    std::vector<std::pair<connection_kind_type, std::size_t>>
    ignore_within() const {return ignore_within_;}

    bool is_ignored_chain(
            const chain_id_type& i, const chain_id_type& j) const noexcept
    {
        return ignored_chain_.is_ignored(i, j);
    }

    static const char* name() noexcept {return "DebyeHuckel";}

    // access to the parameters
    std::vector<real_type>&       charges()       noexcept {return charges_;}
    std::vector<real_type> const& charges() const noexcept {return charges_;}

    //XXX this one is calculated parameter, shouldn't be changed!
    real_type debye_length() const noexcept {return this->debye_length_;}

  private:

    void calc_parameters() noexcept
    {
        constexpr real_type pi =  math::constants<real_type>::pi;
        const real_type kB   = physics::constants<real_type>::kB;
        const real_type e    = physics::constants<real_type>::e;
        const real_type NA   = physics::constants<real_type>::NA;
        const real_type eps0 = physics::constants<real_type>::vacuum_permittivity;
        const real_type epsk = calc_dielectric_water(temperature_, ion_conc_);

        const real_type I = 0.5 * 1000 * ion_conc_;

        this->inv_4_pi_eps0_epsk_ = 1. / (4 * pi * eps0 * epsk);
        this->debye_length_ = std::sqrt(eps0 * epsk * kB * temperature_ /
                                        (2 * NA * e * e * I));
        this->inv_debye_length_ = 1. / this->debye_length_;
        return;
    }

    real_type calc_dielectric_water(
            const real_type T, const real_type c) const noexcept
    {
        return (249.4 - 0.788 * T + 7.2e-4 * T * T) *
               (1. - 2.551 * c + 5.151e-2 * c * c - 6.889e-3 * c * c * c);
    }

  private:

    real_type temperature_;  // [K]
    real_type ion_conc_; // [M]
    real_type inv_4_pi_eps0_epsk_;
    real_type debye_length_;
    real_type inv_debye_length_;
    container_type charges_;

    chain_ignoration_type ignored_chain_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};

template<typename realT, typename ignoreT>
constexpr typename DebyeHuckelPotential<realT, ignoreT>::real_type
DebyeHuckelPotential<realT, ignoreT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_DEBYE_HUCKEL_POTENTIAL */
