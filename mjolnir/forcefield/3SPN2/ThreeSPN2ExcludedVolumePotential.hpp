#ifndef MJOLNIR_FORCEFIELD_3SPN2_EXV_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_3SPN2_EXV_POTENTIAL_HPP
#include <mjolnir/forcefield/3SPN2/ThreeSPN2Common.hpp>
#include <mjolnir/forcefield/global/ParameterList.hpp>
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
    using real_type = realT;

    struct parameter_type
    {
        real_type sigma;
    };

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(1.0);
    }

  public:

    explicit ThreeSPN2ExcludedVolumePotential() noexcept {}
    ~ThreeSPN2ExcludedVolumePotential() = default;

    real_type potential(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma <= r){return 0;}

        const real_type r1s1 = params.sigma / r;
        const real_type r3s3 = r1s1 * r1s1 * r1s1;
        const real_type r6s6 = r3s3 * r3s3;
        return this->epsilon_ * (r6s6 * (r6s6 - real_type(2)) + real_type(1));
    }
    real_type derivative(const real_type r, const parameter_type& params) const noexcept
    {
        if(params.sigma <= r){return 0;}

        const real_type rinv = 1 / r;
        const real_type r1s1 = params.sigma * rinv;
        const real_type r3s3 = r1s1 * r1s1 * r1s1;
        const real_type r6s6 = r3s3 * r3s3;
        return real_type(12) * this->epsilon_ * rinv * r6s6 * (real_type(1) - r6s6);
    }

    template<typename T>
    void initialize(const System<T>& sys) noexcept
    {
        this->update(sys);
        return;
    }

    template<typename T>
    void update(const System<T>&) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        MJOLNIR_LOG_INFO("checking units of parameters...");
        this->epsilon_ = real_type(1); // [kJ/mol]

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
        }
        return;
    }

    template<typename InputIterator>
    real_type max_cutoff(const InputIterator first, const InputIterator last) const noexcept
    {
        static_assert(std::is_same<
                typename std::iterator_traits<InputIterator>::value_type,
                parameter_type>::value, "");

        if(first == last) {return 1;}

        real_type max_sigma = 0;
        for(auto iter = first; iter != last; ++iter)
        {
            const auto& parameter = *iter;
            max_sigma = std::max(max_sigma, parameter.sigma);
        }
        return max_sigma;
    }
    // It returns absolute cutoff length using pair-parameter.
    // `CombinationTable` uses this.
    real_type absolute_cutoff(const parameter_type& params) const noexcept
    {
        return params.sigma;
    }

    static const char* name() noexcept {return "3SPN2ExcludedVolume";}

    real_type epsilon()        const noexcept {return this->epsilon_;}
    real_type cutoff_ratio()   const noexcept {return real_type(1);}
    real_type coef_at_cutoff() const noexcept {return real_type(0);}

  private:

    real_type epsilon_; // 1 in [kJ/mol]
};

// ===========================================================================

template<typename traitsT>
class ThreeSPN2ExcludedVolumeParameterList final
    : public ParameterListBase<traitsT, ThreeSPN2ExcludedVolumePotential<typename traitsT::real_type>>
{
  public:
    using traits_type           = traitsT;
    using real_type             = typename traits_type::real_type;
    using potential_type        = ThreeSPN2ExcludedVolumePotential<real_type>;
    using base_type             = ParameterListBase<traits_type, potential_type>;

    using parameter_type        = parameter_3SPN2::bead_kind;
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
    ThreeSPN2ExcludedVolumeParameterList(const ParameterSet params,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : sigma_P_(params.sigma_P), sigma_S_(params.sigma_S),
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
                this->parameters_.resize(idx+1, parameter_3SPN2::bead_kind::Unknown);
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~ThreeSPN2ExcludedVolumeParameterList() override = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        return pair_parameter_type{(this->sigma_of(this->parameters_[i]) +
                                    this->sigma_of(this->parameters_[j])) * real_type(0.5)};
    }

    real_type max_cutoff_length() const noexcept override
    {
        return this->cutoff_;
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(!unit_converted_)
        {
            MJOLNIR_LOG_INFO("checking units of parameters...");
            // checking unit system and adjust parameters to it

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
        this->update(sys, topol, pot);
        return;
    }

    // nothing to do when system parameters change.
    void update(const system_type& sys, const topology_type& topol,
                const potential_type&) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        using namespace parameter_3SPN2;

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

    // ------------------------------------------------------------------------
    // to check bases has base-pairing interaction.
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        using namespace parameter_3SPN2;
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

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    base_type* clone() const override
    {
        return new ThreeSPN2ExcludedVolumeParameterList(*this);
    }

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "3SPN2ExcludedVolume";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    real_type sigma_of(parameter_3SPN2::bead_kind bk) const noexcept
    {
        using namespace parameter_3SPN2;
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
    real_type epsilon = real_type(1.0); // [kJ/mol]
    real_type sigma_P = real_type(4.5); // [angstrom]
    real_type sigma_S = real_type(6.2);
    real_type sigma_A = real_type(5.4);
    real_type sigma_T = real_type(7.1);
    real_type sigma_G = real_type(4.9);
    real_type sigma_C = real_type(6.4);
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class ThreeSPN2ExcludedVolumeParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ThreeSPN2ExcludedVolumeParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ThreeSPN2ExcludedVolumeParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ThreeSPN2ExcludedVolumeParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif // MJOLNIR_POTENTIAL_GLOBAL_3SPN2_EXV_POTENTIAL_HPP
