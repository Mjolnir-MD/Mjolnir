#ifndef MJOLNIR_FORCEFIELD_3SPN2_EXV_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_3SPN2_EXV_POTENTIAL_HPP
#include <mjolnir/forcefield/3SPN2/ThreeSPN2Common.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstdint>

namespace mjolnir
{

// It calculates a excluded volume that is a part of 3SPN2 DNA model.
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
// Note: an identifier starts with a digit is not allowed in C++ standard.
//       see N3337 2.11 for detail. So `3SPN2BaseBaseInteraction` is not a valid name.
//
template<typename realT>
class ThreeSPN2ExcludedVolumePotential
{
  public:
    using real_type      = realT;
    using parameter_type = real_type; // sigma_ij
    using self_type      = ThreeSPN2ExcludedVolumePotential<real_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(1);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return bead_kind::Unknown;
    }

    static real_type epsilon; // 1 in [kJ/mol]
//     static real_type cutoff_ratio;   // 1 sigma
//     static real_type coef_at_cutoff; // 0.0

  public:

    explicit ThreeSPN2ExcludedVolumePotential(const parameter_type& params) noexcept
        : sigma_(params)
    {}
    ~ThreeSPN2ExcludedVolumePotential() = default;

    real_type potential(const real_type r) const noexcept
    {
        if(sigma_ <= r){return 0;}

        const real_type r1s1   = sigma_ / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        return epsilon * (r6s6 * (r6s6 - real_type(2)) + real_type(1));
    }
    real_type derivative(const real_type r) const noexcept
    {
        if(this->sigma_ <= r){return 0;}

        const real_type rinv   = 1 / r;
        const real_type r1s1   = sigma_ * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        return real_type(12) * epsilon * rinv * r6s6 * (real_type(1) - r6s6);
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "3SPN2ExcludedVolume";}

    real_type sigma()  const noexcept {return this->sigma_;}
    real_type cutoff() const noexcept {return this->sigma_;}

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
            return lhs < rhs;
        }
    };

  private:

    real_type sigma_;
};

// ===========================================================================

template<typename traitsT>
class ThreeSPN2ExcludedVolumeParameter
    : public ParameterListBase<traitsT, ThreeSPN2ExcludedVolumePotential<typename traitsT::real_type>>
{
  public:
    using traits_type           = traitsT;
    using real_type             = typename traits_type::real_type;
    using potential_type        = ThreeSPN2ExcludedVolumePotential<real_type>;
    using base_type             = ParameterListBase<traits_type, potential_type>;

    using parameter_type        = std::pair<real_type, real_type>;
    using pair_parameter_type   = typename potential_type::parameter_type;
    using container_type        = std::vector<parameter_type>;

    // topology stuff
    using system_type           = System<traits_type>;
    using topology_type         = Topology;
    using molecule_id_type      = typename topology_type::molecule_id_type;
    using group_id_type         = typename topology_type::group_id_type;
    using connection_kind_type  = typename topology_type::connection_kind_type;
    using ignore_molecule_type  = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type     = IgnoreGroup   <group_id_type>;
    using exclusion_list_type   = ExclusionList <traits_type>;

  public:

    template<typename ParameterSet>
    ThreeSPN2ExcludedVolumePotential(const ParameterSet params,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : epsilon_(params.epsilon),
          sigma_P_(params.sigma_P), sigma_S_(params.sigma_S),
          sigma_A_(params.sigma_A), sigma_T_(params.sigma_T),
          sigma_G_(params.sigma_G), sigma_C_(params.sigma_C),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // Normally, DNA has two chains and this potential should be applied to
        // both inter-strand and intra-strand particle pairs. So the parameter
        // should be `ignore.molecule = "Nothing"`.
        if(this->exclusion_list_.ignored_molecule_type() != "Nothing"_s)
        {
            MJOLNIR_LOG_WARN("3SPN2 potential requires ignore.molecule = "
                             "\"Nothing\" but you manually set \"",
                             this->exclusion_list_.ignored_molecule_type(),
                             "\". I trust that you know what you are doing.");
        }

        // WCA takes 0 at the sigma. So the cutoff is at the largest sigma.
        this->cutoff_ = std::max(sigma_P_, std::max(sigma_S_, std::max(sigma_A_,
                        std::max(sigma_T_, std::max(sigma_G_, sigma_C_)))));

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
    ~ThreeSPN2ExcludedVolumePotential() = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        return (this->sigma_of(this->parameters_[i]) +
                this->sigma_of(this->parameters_[j])) * 0.5;
    }

    real_type max_cutoff_length() const noexcept override
    {
        return this->cutoff_;
    }

    void initialize(const system_type& sys, const topology_type& topol) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(!unit_converted_)
        {
            MJOLNIR_LOG_INFO("checking units of parameters...");

            // checking unit system and adjust parameters to it
            const auto& energy_unit = physics::constants<real_type>::energy_unit();
            MJOLNIR_LOG_INFO("energy unit is ", energy_unit);
            assert(energy_unit == "kJ/mol" || energy_unit == "kcal/mol");

            if(energy_unit == "kcal/mol")
            {
                MJOLNIR_LOG_INFO("energy unit ([kcal/mol]) differs from the "
                                 "default, [kJ/mol]. converting by multiplying ",
                                 unit::constants<real_type>::J_to_cal());

                // convert from kJ/mol to kcal/mol (/= 4.18)
                this->epsilon_ *= unit::constants<real_type>::J_to_cal();
                ThreeSPN2ExcludedVolumePotential<real_type>::epsilon = epsilon_;
            }

            const auto& length_unit = physics::constants<real_type>::length_unit();
            MJOLNIR_LOG_INFO("length unit is ", length_unit);
            assert(length_unit == "nm" || length_unit == "angstrom");

            if(length_unit == "nm")
            {
                MJOLNIR_LOG_INFO("length unit (nm) differs from the default, "
                                 "[angstrom]. converting by multiplying ",
                                 unit::constants<real_type>::angstrom_to_nm());

                // convert angstrom -> nm (* 0.1)
                this->cutoff_   *= unit::constants<real_type>::angstrom_to_nm();

                this->sigma_P_  *= unit::constants<real_type>::angstrom_to_nm();
                this->sigma_S_  *= unit::constants<real_type>::angstrom_to_nm();
                this->sigma_A_  *= unit::constants<real_type>::angstrom_to_nm();
                this->sigma_T_  *= unit::constants<real_type>::angstrom_to_nm();
                this->sigma_G_  *= unit::constants<real_type>::angstrom_to_nm();
                this->sigma_C_  *= unit::constants<real_type>::angstrom_to_nm();
            }
            unit_converted_ = true;
        }
        this->update(sys, topol);
        return;
    }

    // nothing to do when system parameters change.
    void update(const system_type& sys, const topology_type& topol) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // make exclusion list based on the topology
        exclusion_list_.make(sys, topol);

        // --------------------------------------------------------------------
        // list up beads that are within 3 nucleotides
        this->within_3_nucl_.reserve(this->participants_.size() * 2);

        const auto is_base = [this](const std::size_t i) noexcept -> bool {
            const auto bk = this->parameters_.at(i);
            return bk == bead_kind::Phosphate || bk == bead_kind::Sugar;
        };
        std::ptrdiff_t idx = 0;
        for(const auto i : this->participants_)
        {
            const auto first = idx;
            if(!is_base(i))
            {
                // if it is not a base, we don't need the list.
                // fill it by a empty range and skip this.
                this->within_3_nucl_ranges_.emplace_back(first, first);
                continue;
            }
            auto within3 = topol.list_adjacent_within(i, 3, "next_nucleotide");

            // remove non-base stuff
            within3.erase(std::remove_if(within3.begin(), within3.end(), is_base),
                          within3.end());
            // make it unique
            std::sort(within3.begin(), within3.end());
            within3.erase(std::unique(within3.begin(), within3.end()),
                          within3.end());

            // append it to the internal cache
            for(const auto j : within3)
            {
                this->within_3_nucl_.push_back(j);
                ++idx;
            }
            // and keep the range corresponds to i to the container of ranges.
            this->within_3_nucl_ranges_.emplace_back(first, idx);
        }
        return;
    }

    // -----------------------------------------------------------------------
    // for spatial partitions

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

    // ------------------------------------------------------------------------
    // to check bases has base-pairing interaction.
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        if(j <= i || exclusion_list_.is_excluded(i, j))
        {
            return false;
        }

        const auto& i_kind = this->parameters_[i];
        const auto& j_kind = this->parameters_[j];
        assert(i_kind != bead_kind::Unknown);
        assert(j_kind != bead_kind::Unknown);

        // if not a base, they always has interaction.
        if(i_kind == bead_kind::Phosphate || i_kind == bead_kind::Sugar ||
           j_kind == bead_kind::Phosphate || j_kind == bead_kind::Sugar)
        {
            return true;
        }

        // If bases are not separated by 3 nucleotides, they cannot form
        // base-pairing. In that case, excluded volume will be applied.
        const auto within3 = this->within_3_nucleotides(i);
        if(std::binary_search(within3.begin(), within3.end(), j))
        {
            return true;
        }

        // if the bases can form a base-pairing, then it does not have interaction.
        switch(i_kind)
        {
            case bead_kind::BaseA: {return j_kind != bead_kind::BaseT;}
            case bead_kind::BaseT: {return j_kind != bead_kind::BaseA;}
            case bead_kind::BaseG: {return j_kind != bead_kind::BaseC;}
            case bead_kind::BaseC: {return j_kind != bead_kind::BaseG;}
            default: {assert(false);}
        }
    }

    exclusion_list_type const& exclusion_list() const noexcept {return exclusion_list_;}

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "3SPN2ExcludedVolume";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    real_type sigma_of(bead_kind bk) const noexcept
    {
        switch(bk)
        {
            case bead_kind::Phosphate: {return this->sigma_P_;}
            case bead_kind::Sugar:     {return this->sigma_S_;}
            case bead_kind::BaseA:     {return this->sigma_A_;}
            case bead_kind::BaseT:     {return this->sigma_T_;}
            case bead_kind::BaseC:     {return this->sigma_C_;}
            case bead_kind::BaseG:     {return this->sigma_G_;}
            default:{assert(false);}
        }
    }

  private:

    range<typename std::vector<std::size_t>::const_iterator>
    within_3_nucleotides(const std::size_t i) const noexcept
    {
        return range<typename std::vector<std::size_t>::const_iterator>{
            this->within_3_nucl_.begin() + this->within_3_nucl_ranges_[i].first,
            this->within_3_nucl_.begin() + this->within_3_nucl_ranges_[i].second
        };
    }

  private:

    bool unit_converted_ = false;
    real_type epsilon_; // [kJ/mol]
    real_type sigma_P_; // [angstrom]
    real_type sigma_S_;
    real_type sigma_A_;
    real_type sigma_T_;
    real_type sigma_G_;
    real_type sigma_C_;
    real_type cutoff_; // the maximum sigma_ij

    container_type           parameters_;
    std::vector<std::size_t> participants_;
    exclusion_list_type      exclusion_list_;

    // list of bases within 3 nucleotides.
    //
    // If bases are not separated by 3 nucleotides, they cannot form base pair.
    // In that case, excluded volume will be applied.
    std::vector<std::size_t> within_3_nucl_;
    std::vector<std::pair<std::ptrdiff_t, std::ptrdiff_t>> within_3_nucl_ranges_;
};

template<typename realT>
struct ThreeSPN2ExcludedVolumePotentialParameter
{
    using real_type = realT;
    real_type epsilon = 1.0; // [kJ/mol]
    real_type sigma_P = 4.5; // [angstrom]
    real_type sigma_S = 6.2;
    real_type sigma_A = 5.4;
    real_type sigma_T = 7.1;
    real_type sigma_G = 4.9;
    real_type sigma_C = 6.4;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ThreeSPN2ExcludedVolumePotential<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif // MJOLNIR_POTENTIAL_GLOBAL_3SPN2_EXV_POTENTIAL_HPP
