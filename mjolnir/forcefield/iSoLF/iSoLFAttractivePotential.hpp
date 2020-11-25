#ifndef MJOLNIR_POTENTIAL_GLOBAL_ISOLF_ATTRACTIVE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_ISOLF_ATTRACTIVE_POTENTIAL_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <cmath>

namespace mjolnir
{

// attractive part of the iSoLF potential for coarse-grained lipids developed by
// - Diego Ugarte La Torre and Shoji Takada (2020) J. Chem. Phys 153, 205101
//   https://doi.org/10.1063/5.0026342
//
// This class contains sigmas, epsilons, and omegas of the particles and
// calculates energy and derivative of the potential function.
template<typename traitsT>
class iSoLFAttractivePotential
{
  public:
    using traits_type          = traitsT;
    using real_type            = typename traits_type::real_type;
    using system_type          = System<traits_type>;
    using parameter_type       = std::tuple<real_type, real_type, real_type>;
    using container_type       = std::vector<parameter_type>;

    // cache sigma_ij, epsilon_ij, omega_ij, 1 / 2omega_ij.
    using pair_parameter_type  = std::tuple<real_type, real_type, real_type, real_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0), real_type(0)};
    }

  public:

    iSoLFAttractivePotential(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
    : exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
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
    }
    ~iSoLFAttractivePotential() = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        const auto sgm1 = std::get<0>(parameters_[i]);
        const auto eps1 = std::get<1>(parameters_[i]);
        const auto omg1 = std::get<2>(parameters_[i]);

        const auto sgm2 = std::get<0>(parameters_[j]);
        const auto eps2 = std::get<1>(parameters_[j]);
        const auto omg2 = std::get<2>(parameters_[j]);

        return std::make_tuple((sgm1 + sgm2) / 2,
                              ((eps1 == eps2) ? eps1 : std::sqrt(eps1 * eps2)),
                               (omg1 + omg2) / 2,
                               real_type(1.0) / (omg1 + omg2));
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

    real_type potential(const real_type r, const pair_parameter_type& p) const noexcept
    {
        constexpr real_type rc = 1.12246204831;
        constexpr real_type pi = math::constants<real_type>::pi();

        const real_type sigma   = std::get<0>(p);
        const real_type epsilon = std::get<1>(p);
        const real_type omega   = std::get<2>(p);

        const real_type r_sigma_rc = r - sigma * rc; // r - sqrt[6]{2} sigma

        if     (r_sigma_rc <= 0)   {return -epsilon;}
        else if(omega < r_sigma_rc){return 0;}

        const real_type romega = std::get<3>(p); // 1 / 2omega_ij
        const real_type cosine = std::cos(pi * romega * r_sigma_rc);

        return -epsilon * cosine * cosine;
    }
    real_type derivative(const real_type r, const pair_parameter_type& p) const noexcept
    {
        constexpr real_type rc = 1.12246204831;
        constexpr real_type pi = math::constants<real_type>::pi();

        const real_type sigma   = std::get<0>(p);
        const real_type epsilon = std::get<1>(p);
        const real_type omega   = std::get<2>(p);

        const real_type r_sigma_rc = r - sigma * rc; // r - sqrt[6]{2} sigma

        if (r_sigma_rc <= 0 || omega < r_sigma_rc) {return 0;}

        const real_type romega = std::get<3>(p); // 1 / 2omega_ij
        const real_type sine   = std::sin(2 * pi * romega * r_sigma_rc);

        return epsilon * pi * romega * sine;
    }

    real_type cutoff_ratio()   const noexcept {return std::numeric_limits<real_type>::infinity();}
    real_type coef_at_cutoff() const noexcept {return 0.0;}

    real_type max_cutoff_length() const noexcept
    {
        constexpr real_type rc = 1.12246204831; // sqrt[6]{2}

        if(this->parameters_.empty())
        {
            return 0.0;
        }

        const real_type max_sigma = std::get<0>(*std::max_element(
            this->parameters_.cbegin(), this->parameters_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return std::get<0>(lhs) < std::get<0>(rhs);
            }));
        const real_type max_omega = std::get<2>(*std::max_element(
            this->parameters_.cbegin(), this->parameters_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return std::get<2>(lhs) < std::get<2>(rhs);
            }));

        return max_sigma * rc + max_omega;
    }

    void initialize(const system_type& sys, const topology_type& topol) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol);
        return;
    }

    void update(const system_type& sys, const topology_type& topol) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

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
        return make_range(participants_.begin() + participant_idx + 1, participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "iSoLF";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

  private:

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
extern template class iSoLFAttractivePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class iSoLFAttractivePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class iSoLFAttractivePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class iSoLFAttractivePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_ISOLF_ATTRACTIVE_POTENTIAL */
