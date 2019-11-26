#ifndef MJOLNIR_POTENTIAL_GLOBAL_UNIFORM_LENNARD_JONES_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_UNIFORM_LENNARD_JONES_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/empty.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace mjolnir
{

// Well-known Lennard-Jones interaction with uniform parameters.
// This class contains a sigma and an epsilon that are the same among all the
// particles.
template<typename traitsT>
class UniformLennardJonesPotential
{
  public:
    using traits_type         = traitsT;
    using real_type           = typename traits_type::real_type;
    using parameter_type      = empty_t; // no particle-specific parameter
    using container_type      = empty_t; // no parameter, so no container there.
    using pair_parameter_type = empty_t; // no particle-pair-specific parameter

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.5);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{};
    }

  public:

    UniformLennardJonesPotential(
        const real_type sgm, const real_type eps, const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
    : cutoff_ratio_(cutoff_ratio),
      coef_at_cutoff_(std::pow(1/cutoff_ratio, 12) - std::pow(1/cutoff_ratio, 6)),
      sigma_(sgm), epsilon_(eps), r_cut_(sgm * cutoff_ratio),
      exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            this->participants_.push_back(idxp.first);
        }
    }
    ~UniformLennardJonesPotential() = default;

    pair_parameter_type prepare_params(std::size_t, std::size_t) const noexcept
    {
        return pair_parameter_type{}; // no pre-calculated parameter
    }

    // forwarding functions for clarity...
    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const noexcept
    {
        return this->potential(r, this->prepare_params(i, j));
    }
    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        return this->derivative(r, this->prepare_params(i, j));
    }

    real_type
    potential(const real_type r, const pair_parameter_type&) const noexcept
    {
        if(r_cut_ < r){return 0;}

        const real_type r1s1   = sigma_ / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 4 * epsilon_ * (r12s12 - r6s6 - this->coef_at_cutoff_);
    }

    real_type
    derivative(const real_type r, const pair_parameter_type&) const noexcept
    {
        if(r_cut_ < r){return 0;}

        const real_type rinv   = 1 / r;
        const real_type r1s1   = sigma_ * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 24 * epsilon_ * (r6s6 - 2 * r12s12) * rinv;
    }

    real_type cutoff_ratio()   const noexcept {return this->cutoff_ratio_;}
    real_type coef_at_cutoff() const noexcept {return this->coef_at_cutoff_;}

    real_type max_cutoff_length() const noexcept
    {
        return sigma_ * cutoff_ratio_;
    }

    void initialize(const System<traits_type>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // if no participants are given, consider all the particles are related.
        if(this->participants_.empty())
        {
            MJOLNIR_LOG_WARN("UniformLennardJonesPotential does not have any participants.");
            MJOLNIR_LOG_WARN("All the particles are considered to be participants.");
            MJOLNIR_LOG_WARN("To disable this potential, "
                             "simply remove the part from the input file.");

            this->participants_.resize(sys.size());
            std::iota(this->participants_.begin(), this->participants_.end(), 0u);
        }

        this->update(sys);
        return;
    }

    void update(const System<traits_type>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys);
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

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}

    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept
    {
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "LennardJones";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    real_type& sigma()         noexcept {return sigma_;}
    real_type  sigma()   const noexcept {return sigma_;}
    real_type& epsilon()       noexcept {return epsilon_;}
    real_type  epsilon() const noexcept {return epsilon_;}

  private:

    real_type cutoff_ratio_, coef_at_cutoff_;
    real_type sigma_, epsilon_, r_cut_;
    std::vector<std::size_t> participants_;

    exclusion_list_type  exclusion_list_;
};
} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class UniformLennardJonesPotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class UniformLennardJonesPotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class UniformLennardJonesPotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class UniformLennardJonesPotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
