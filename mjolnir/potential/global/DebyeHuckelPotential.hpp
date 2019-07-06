#ifndef MJOLNIR_POTENTIAL_GLOBAL_DEBYE_HUCKEL_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_DEBYE_HUCKEL_POTENTIAL_HPP
#include <mjolnir/potential/global/IgnoreMolecule.hpp>
#include <mjolnir/potential/global/IgnoreGroup.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <numeric>
#include <cmath>
#include <cassert>

namespace mjolnir
{

// Debye-Huckel type electrostatic interaction.
// This class contains charges and other parameters and calculates energy and
// derivative of the potential function.
// This class is an implementation of the electrostatic term used in the 3SPN
// series of the coarse-grained DNA models (Knotts et al., (2007), Sambriski and
// de Pablo, (2009), Hinckley et al., (2013), and Freeman et al., (2014))
// It is same as the default electrostatic term in CafeMol (Kenzaki et al. 2011)
template<typename realT>
class DebyeHuckelPotential
{
  public:
    using real_type            = realT;
    using parameter_type       = real_type;
    using container_type       = std::vector<parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;

    // r_cutoff = cutoff_ratio * debye_length
    static constexpr real_type cutoff_ratio = 5.5;

    static constexpr parameter_type default_parameter() noexcept
    {
        return real_type(0);
    }

  public:

    DebyeHuckelPotential(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : temperature_(300.0), ion_conc_(0.1),
          ignore_molecule_(std::move(ignore_mol)),
          ignore_group_   (std::move(ignore_grp)),
          ignore_within_  (exclusions.begin(), exclusions.end())
    {
        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }

        // XXX should be updated before use because T and ion conc are default!
        this->calc_parameters();
    }
    ~DebyeHuckelPotential() = default;

    parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return this->inv_4_pi_eps0_epsk_ * this->parameters_[i] * this->parameters_[j];
    }

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type   r) const noexcept
    {
        return this->potential(r, this->prepare_params(i, j));
    }
    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type   r) const noexcept
    {
        return this->derivative(r, this->prepare_params(i, j));
    }

    real_type potential(const real_type r, const parameter_type& p) const noexcept
    {
        if(this->max_cutoff_length() <= r) {return 0.0;}
        return p * std::exp(-r * this->inv_debye_length_) / r;
    }
    real_type derivative(const real_type r, const parameter_type& p) const noexcept
    {
        if(this->max_cutoff_length() <= r) {return 0.0;}
        return -p * (debye_length_ + r) * this->inv_debye_length_ / (r * r) *
               std::exp(-r * this->inv_debye_length_);
    }

    real_type max_cutoff_length() const noexcept
    {
        return debye_length_ * cutoff_ratio;
    }

    template<typename traitsT>
    void initialize(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(!sys.has_attribute("temperature"))
        {
            MJOLNIR_LOG_ERROR("DebyeHuckel requires `temperature` attribute");
        }
        if(!sys.has_attribute("ionic_strength"))
        {
            MJOLNIR_LOG_ERROR("DebyeHuckel requires `ionic_strength` attribute");
        }
        this->temperature_ = sys.attribute("temperature");
        this->ion_conc_    = sys.attribute("ionic_strength");
        this->calc_parameters();
        return;
    }

    // for temperature/ionic concentration changes...
    template<typename traitsT>
    void update(const System<traitsT>& sys) noexcept
    {
        assert(sys.has_attribute("temperature"));
        assert(sys.has_attribute("ionic_strength"));

        this->temperature_ = sys.attribute("temperature");
        this->ion_conc_    = sys.attribute("ionic_strength");
        this->calc_parameters();
        return;
    }

    // e.g. {"bond", 3} means ignore particles connected within 3 "bond"s
    std::vector<std::pair<connection_kind_type, std::size_t>>
    ignore_within() const {return ignore_within_;}

    bool is_ignored_molecule(
            const molecule_id_type& i, const molecule_id_type& j) const noexcept
    {
        return ignore_molecule_.is_ignored(i, j);
    }
    bool is_ignored_group(
            const group_id_type& i, const group_id_type& j) const noexcept
    {
        return ignore_group_.is_ignored(i, j);
    }

    static const char* name() noexcept {return "DebyeHuckel";}

    // access to the parameters
    std::vector<real_type>&       charges()       noexcept {return parameters_;}
    std::vector<real_type> const& charges() const noexcept {return parameters_;}

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}

    //XXX this one is calculated parameter, shouldn't be changed!
    real_type debye_length() const noexcept {return this->debye_length_;}

  private:

    void calc_parameters() noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        using math_const =    math::constants<real_type>;
        using phys_const = physics::constants<real_type>;

        constexpr real_type pi   = math_const::pi;
        const     real_type kB   = phys_const::kB();
        const     real_type NA   = phys_const::NA();
        const     real_type eps0 = phys_const::eps0();

        const     real_type epsk = calc_dielectric_water(temperature_, ion_conc_);
        const     real_type T    = this->temperature_;

        MJOLNIR_LOG_INFO("kB               = ", kB);
        MJOLNIR_LOG_INFO("T                = ", T);
        MJOLNIR_LOG_INFO("NA               = ", NA);
        MJOLNIR_LOG_INFO("epsilon_0        = ", eps0);
        MJOLNIR_LOG_INFO("epsilon_k        = ", epsk);

        this->inv_4_pi_eps0_epsk_ = 1.0 / (4 * pi * eps0 * epsk);

        // convert [M] (mol/L) or [mM] -> [mol/nm^3] or [mol/A^3]
        const real_type I = 0.5 * ion_conc_ / phys_const::L_to_volume();

        this->debye_length_ = std::sqrt((eps0 * epsk * kB * T) / (2 * NA * I));
        this->inv_debye_length_ = 1. / this->debye_length_;

        MJOLNIR_LOG_INFO("debye length      = ", debye_length_);
        MJOLNIR_LOG_INFO("1 / debye length  = ", inv_debye_length_);
        MJOLNIR_LOG_INFO("1 / 4pi eps0 epsk = ", inv_4_pi_eps0_epsk_);

        return;
    }

    real_type calc_dielectric_water(
            const real_type T, const real_type c) const noexcept
    {
        return (249.4 - 0.788 * T + 7.2e-4 * T * T) *
               (1. - 2.551e-1 * c + 5.151e-2 * c * c - 6.889e-3 * c * c * c);
    }

  private:

    real_type temperature_;  // [K]
    real_type ion_conc_; // [M]
    real_type inv_4_pi_eps0_epsk_;
    real_type debye_length_;
    real_type inv_debye_length_;

    container_type parameters_;
    std::vector<std::size_t> participants_;

    ignore_molecule_type ignore_molecule_;
    ignore_group_type    ignore_group_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};

template<typename realT>
constexpr typename DebyeHuckelPotential<realT>::real_type
DebyeHuckelPotential<realT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_DEBYE_HUCKEL_POTENTIAL */
