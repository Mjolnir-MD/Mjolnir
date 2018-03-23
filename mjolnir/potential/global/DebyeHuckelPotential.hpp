#ifndef MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#define MJOLNIR_DEBYE_HUCKEL_POTENTIAL
#include <mjolnir/core/constants.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/potential/global/GroupIgnoration.hpp>
#include <vector>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, template<typename GID> class GroupIgnoration>
class DebyeHuckelPotential
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef real_type parameter_type;
    typedef std::vector<parameter_type> container_type;

    // topology stuff
    typedef StructureTopology topology_type;
    typedef typename topology_type::group_id_type        group_id_type;
    typedef typename topology_type::connection_kind_type connection_kind_type;
    typedef GroupIgnoration<group_id_type> group_ignoration_type;

    // r_cutoff = cutoff_ratio * debye_length
    constexpr static real_type cutoff_ratio = 5.5;

  public:

    DebyeHuckelPotential() = default;
    DebyeHuckelPotential(const container_type& charges)
      : charges_(charges), temperature_(300.0), ion_conc_(0.1)
    {
        // XXX should be updated before use because T and ion conc are default!
        this->calc_parameters();
    }
    DebyeHuckelPotential(container_type&& charges)
      : charges_(std::move(charges)), temperature_(300.0), ion_conc_(0.1)
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
    void update(const System<traitsT>& sys) noexcept
    {
        this->temperature_ = sys.attribute("temperature");
        this->ion_conc_    = sys.attribute("ionic_strength");
        this->calc_parameters();
        return;
    }

    // e.g. {"bond", 3} means ignore particles connected within 3 "bond"s
    std::size_t  ignored_bonds()    const noexcept {return ignored_bonds_;}
    std::size_t& ignored_bonds()          noexcept {return ignored_bonds_;}
    std::size_t  ignored_contacts() const noexcept {return ignored_contacts_;}
    std::size_t& ignored_contacts()       noexcept {return ignored_contacts_;}

    bool is_ignored_group(
            const group_id_type& i, const group_id_type& j) const noexcept
    {
        return ignored_group_.is_ignored(i, j);
    }

    std::string name() const noexcept {return "DebyeHuckel";}

    // access to the parameters
    std::vector<real_type>&       charges()       noexcept {return charges_;}
    std::vector<real_type> const& charges() const noexcept {return charges_;}

    //XXX this one is calculated parameter, shouldn't be changed!
    real_type debye_length() const noexcept {return this->debye_length_;}

  private:

    void calc_parameters() noexcept
    {
        const real_type kB   = physics<real_type>::kB;
        const real_type e    = physics<real_type>::e;
        const real_type NA   = physics<real_type>::NA;
        const real_type eps0 = physics<real_type>::vacuum_permittivity;
        const real_type epsk = calc_dielectric_water(temperature_, ion_conc_);
        const real_type pi   = constants<real_type>::pi;

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

    group_ignoration_type ignored_group_;
    std::size_t ignored_bonds_;
    std::size_t ignored_contacts_;
};

template<typename traitsT, template<typename> class ignoreT>
constexpr typename DebyeHuckelPotential<traitsT, ignoreT>::real_type
DebyeHuckelPotential<traitsT, ignoreT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_DEBYE_HUCKEL_POTENTIAL */
