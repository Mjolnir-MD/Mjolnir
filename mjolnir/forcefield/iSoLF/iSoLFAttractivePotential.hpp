#ifndef MJOLNIR_POTENTIAL_GLOBAL_ISOLF_ATTRACTIVE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_ISOLF_ATTRACTIVE_POTENTIAL_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
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
template<typename realT>
class iSoLFAttractivePotential
{
  public:
    using real_type      = realT; // sigma, epsilon, omega
    using parameter_type = std::tuple<real_type, real_type, real_type>;
    using self_type      = iSoLFAttractivePotential<real_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return std::numeric_limits<real_type>::infinity();;
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0), real_type(0)};
    }

    static void set_cutoff_ratio(const real_type)
    {
        return;
    }

  public:

    iSoLFAttractivePotential(const parameter_type& params) noexcept
        : sigma_(std::get<0>(params)), epsilon_(std::get<1>(params)),
          omega_(std::get<2>(params)), romega_(real_type(1) / std::get<2>(params))
    {}
    ~iSoLFAttractivePotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        constexpr real_type rc = 1.12246204831;
        constexpr real_type pi = math::constants<real_type>::pi();

        const real_type r_sigma_rc = r - sigma_ * rc; // r - sqrt[6]{2} sigma

        if     (r_sigma_rc <= 0)    {return -epsilon_;}
        else if(omega_ < r_sigma_rc){return 0;}

        const real_type cosine = std::cos(pi * romega_ * r_sigma_rc);

        return -epsilon_ * cosine * cosine;
    }
    real_type derivative(const real_type r) const noexcept
    {
        constexpr real_type rc = 1.12246204831;
        constexpr real_type pi = math::constants<real_type>::pi();

        const real_type r_sigma_rc = r - sigma_ * rc; // r - sqrt[6]{2} sigma

        if (r_sigma_rc <= 0 || omega_ < r_sigma_rc) {return 0;}

        const real_type sine = std::sin(2 * pi * romega_ * r_sigma_rc);

        return epsilon_ * pi * romega_ * sine;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "iSoLFAttractive";}

    real_type sigma()   const noexcept {return this->sigma_;}
    real_type epsilon() const noexcept {return this->epsilon_;}
    real_type omega()   const noexcept {return this->omega_;}

    real_type cutoff()  const noexcept
    {
        constexpr real_type rc = 1.12246204831;
        return sigma_ * rc + omega_;
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
            constexpr real_type rc = 1.12246204831;
            return std::get<0>(lhs) * rc + std::get<1>(lhs) <
                   std::get<0>(rhs) * rc + std::get<1>(rhs) ;
        }
    };

  private:

    real_type sigma_;
    real_type epsilon_;
    real_type omega_;
    real_type romega_;
};

template<typename traitsT>
class iSoLFAttractiveParameterList final
    : public ParameterListBase<traitsT, iSoLFAttractivePotential<typename traitsT::real_type>>
{
  public:
    using traits_type          = traitsT;
    using real_type            = typename traits_type::real_type;
    using potential_type       = iSoLFAttractivePotential<real_type>;
    using base_type            = ParameterListBase<traits_type, potential_type>;

    using parameter_type       = std::tuple<real_type, real_type, real_type>; // sigma, epsilon, omega
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
    using exclusion_list_type  = ExclusionList <traits_type>;

    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0), real_type(0)};
    }

  public:

    iSoLFAttractiveParameterList(
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
    ~iSoLFAttractiveParameterList() = default;

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
                               (omg1 + omg2) / 2);
    }
    real_type max_cutoff_length() const noexcept override
    {
        return this->max_cutoff_length_;
    }

    real_type cutoff_ratio()   const noexcept {return std::numeric_limits<real_type>::infinity();}
    real_type coef_at_cutoff() const noexcept {return 0.0;}

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
        constexpr real_type rc = 1.12246204831; // sqrt[6]{2}

        if(this->parameters_.empty())
        {
            this->max_cutoff_length_ = 1.0;
        }
        else
        {
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

            this->max_cutoff_length_ = max_sigma * rc + max_omega;
        }
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
        return make_range(participants_.begin() + participant_idx + 1, participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    // for testing
    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "iSoLFAttractive";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    base_type* clone() const override
    {
        return new iSoLFAttractiveParameterList(*this);
    }

  private:

    real_type max_cutoff_length_;
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
extern template class iSoLFAttractiveParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class iSoLFAttractiveParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class iSoLFAttractiveParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class iSoLFAttractiveParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_ISOLF_ATTRACTIVE_POTENTIAL */
