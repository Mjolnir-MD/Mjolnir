#ifndef MJOLNIR_POTENTIAL_GLOBAL_DEBYE_HUCKEL_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_DEBYE_HUCKEL_POTENTIAL_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <numeric>
#include <cmath>
#include <cassert>

namespace mjolnir
{

template<typename realT>
class DebyeHuckelPotential
{
  public:
    using real_type      = realT;
    using parameter_type = real_type; // qi qj / 4pi epsilon
    using self_type      = DebyeHuckelPotential<real_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(5.5);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0), real_type(0)};
    }

    static void set_cutoff_ratio(const real_type ratio)
    {
        if(self_type::cutoff_ratio < ratio)
        {
            self_type::cutoff_ratio   = ratio;
            self_type::coef_at_cutoff = std::exp(-self_type::cutoff_ratio) /
                    (self_type::debye_length * self_type::cutoff_ratio);
        }
        return;
    }

    static real_type cutoff_ratio;
    static real_type coef_at_cutoff;
    static real_type debye_length;
    static real_type inv_debye_length;

  public:

    explicit DebyeHuckelPotential(const parameter_type& params) noexcept
        : k_(params)
    {}
    ~DebyeHuckelPotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        if(self_type::cutoff_ratio <= r * inv_debye_length)
        {
            return 0.0;
        }
        return k_ * (std::exp(-r * inv_debye_length) / r - coef_at_cutoff);
    }
    real_type derivative(const real_type r) const noexcept
    {
        if(self_type::cutoff_ratio <= r * inv_debye_length)
        {
            return 0.0;
        }
        return -k_ * (debye_length + r) * inv_debye_length / (r * r) *
               std::exp(-r * inv_debye_length);
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "DebyeHuckel";}

    real_type k() const noexcept {return this->k_;}

    real_type cutoff()  const noexcept
    {
        return debye_length * self_type::cutoff_ratio;
    }

  public:

    // To culculate cutoff distance, we need to find the maximum sigma in the
    // existing parameters. But the list of parameters will be given in a variety
    // of ways, like Lorentz-Bertherot rule, combination table, or another way
    // of combination rules.
    //     To find the maximum parameter, we need to provide a way to compare
    // parameters. But the way depends on the functional form of a potential.
    // So this comparator should be defined in a Potential class.
    struct parameter_comparator
    {
        constexpr bool
        operator()(const parameter_type& lhs, const parameter_type& rhs) const noexcept
        {
            return lhs < rhs; // anything is okay here. this is not used
            // because cutoff length depends only on the debye length, not charges.
        }
    };

  private:

    real_type k_; // qiqj / 4pi epsilon
};

template<typename realT>
typename DebyeHuckelPotential<realT>::real_type DebyeHuckelPotential<realT>::cutoff_ratio  = 0.0;
template<typename realT>
typename DebyeHuckelPotential<realT>::real_type DebyeHuckelPotential<realT>::coef_at_cutoff = 0.0;

template<typename realT>
typename DebyeHuckelPotential<realT>::real_type DebyeHuckelPotential<realT>::debye_length = 0.0;
template<typename realT>
typename DebyeHuckelPotential<realT>::real_type DebyeHuckelPotential<realT>::inv_debye_length = 0.0;

// Debye-Huckel type electrostatic interaction.
// This class contains charges and other parameters and calculates pair parameter
// from those information.
// This class is an implementation of the electrostatic term used in the 3SPN
// series of the coarse-grained DNA models (Knotts et al., (2007), Sambriski and
// de Pablo, (2009), Hinckley et al., (2013), and Freeman et al., (2014))
// It is same as the default electrostatic term in CafeMol (Kenzaki et al. 2011)
template<typename traitsT>
class DebyeHuckelParameterList final
    : public ParameterListBase<traitsT, DebyeHuckelPotential<typename traitsT::real_type>>
{
  public:
    using traits_type          = traitsT;
    using real_type            = typename traits_type::real_type;
    using potential_type       = DebyeHuckelPotential<real_type>;
    using base_type            = ParameterListBase<traits_type, potential_type>;

    // `pair_parameter_type` is a parameter for an interacting pair.
    // Although it is the same type as `parameter_type` in this potential,
    // it can be different from normal parameter for each particle.
    // This enables NeighborList to cache a value to calculate forces between
    // the particles, e.g. by having qi * qj for pair of particles i and j.

    using parameter_type       = real_type; // q_i
    using pair_parameter_type  = typename potential_type::parameter_type;
    using container_type       = std::vector<parameter_type>;

    // topology stuff
    using system_type          = System<traits_type>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList<traits_type>;

  public:

    DebyeHuckelParameterList(const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
    : cutoff_ratio_(cutoff_ratio),
        // XXX should be updated in the `initialize(sys)` method
      temperature_ (std::numeric_limits<real_type>::quiet_NaN()),
      ion_strength_(std::numeric_limits<real_type>::quiet_NaN()),
      exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, real_type(0));
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~DebyeHuckelParameterList() = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        return this->inv_4_pi_eps0_epsk_ * this->parameters_[i] * this->parameters_[j];
    }

    real_type max_cutoff_length() const noexcept override
    {
        return this->debye_length_ * this->cutoff_ratio_;
    }

    void initialize(const system_type& sys, const topology_type& topol) noexcept override
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
        this->update(sys, topol); // calc parameters
        return;
    }

    // for temperature/ionic concentration changes...
    void update(const system_type& sys, const topology_type& topol) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        assert(sys.has_attribute("temperature"));
        assert(sys.has_attribute("ionic_strength"));

        this->temperature_  = sys.attribute("temperature");
        this->ion_strength_ = sys.attribute("ionic_strength");

        MJOLNIR_LOG_INFO("temperature    = ", this->temperature_);
        MJOLNIR_LOG_INFO("ionic strength = ", this->ion_strength_);

        this->calc_parameters();

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys, topol);
        return;
    }

    // -----------------------------------------------------------------------
    // for spatial partitions
    //
    // Here, the default implementation uses Newton's 3rd law to reduce
    // calculation. For an interacting pair (i, j), forces applied to i and j
    // are equal in magnitude and opposite in direction. So, if a pair (i, j) is
    // listed, (j, i) is not needed.
    //     See implementation of VerletList, CellList and GlobalPairInteraction
    // for more details about the usage of these functions.

    std::vector<std::size_t> const& participants() const noexcept override
    {
        return participants_;
    }
    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept override
    {
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_; // for testing
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "DebyeHuckel";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters.
    std::vector<real_type>&       charges()       noexcept {return parameters_;}
    std::vector<real_type> const& charges() const noexcept {return parameters_;}

    real_type cutoff_ratio() const noexcept {return this->cutoff_ratio_;}
    real_type debye_length() const noexcept {return this->debye_length_;}

    base_type* clone() const override
    {
        return new DebyeHuckelParameterList(*this);
    }

  private:

    void calc_parameters() noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        using math_const =    math::constants<real_type>;
        using phys_const = physics::constants<real_type>;

        constexpr real_type pi   = math_const::pi();
        const     real_type kB   = phys_const::kB();
        const     real_type NA   = phys_const::NA();
        const     real_type eps0 = phys_const::eps0();

        const     real_type epsk = calc_dielectric_water(temperature_, ion_strength_);
        const     real_type T    = this->temperature_;

        MJOLNIR_LOG_INFO("kB               = ", kB);
        MJOLNIR_LOG_INFO("T                = ", T);
        MJOLNIR_LOG_INFO("NA               = ", NA);
        MJOLNIR_LOG_INFO("epsilon_0        = ", eps0);
        MJOLNIR_LOG_INFO("epsilon_k        = ", epsk);

        this->inv_4_pi_eps0_epsk_ = 1.0 / (4 * pi * eps0 * epsk);

        // convert [M] (mol/L) or [mM] -> [mol/nm^3] or [mol/A^3]
        const real_type I = ion_strength_ / phys_const::L_to_volume();

        this->debye_length_ = std::sqrt((eps0 * epsk * kB * T) / (2 * NA * I));
        this->inv_debye_length_ = 1. / this->debye_length_;

        MJOLNIR_LOG_INFO("debye length      = ", debye_length_);
        MJOLNIR_LOG_INFO("1 / debye length  = ", inv_debye_length_);
        MJOLNIR_LOG_INFO("1 / 4pi eps0 epsk = ", inv_4_pi_eps0_epsk_);

        // set parameter of potential_type
        potential_type::debye_length     = this->debye_length_;
        potential_type::inv_debye_length = this->inv_debye_length_;
        potential_type::set_cutoff_ratio(cutoff_ratio_);

        return;
    }

    real_type calc_dielectric_water(
            const real_type T, const real_type c) const noexcept
    {
        return (249.4 - 0.788 * T + 7.2e-4 * T * T) *
               (1. - 2.551e-1 * c + 5.151e-2 * c * c - 6.889e-3 * c * c * c);
    }

  private:

    real_type cutoff_ratio_; // relative to the debye length
    real_type temperature_;  // [K]
    real_type ion_strength_; // [M]
    real_type inv_4_pi_eps0_epsk_;
    real_type debye_length_;
    real_type inv_debye_length_;

    container_type parameters_;
    std::vector<std::size_t> participants_;

    exclusion_list_type  exclusion_list_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class DebyeHuckelParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class DebyeHuckelParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class DebyeHuckelParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class DebyeHuckelParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_DEBYE_HUCKEL_POTENTIAL */
