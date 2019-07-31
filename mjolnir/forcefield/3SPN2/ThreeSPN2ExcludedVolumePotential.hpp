#ifndef MJOLNIR_POTENTIAL_GLOBAL_3SPN2_EXV_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_3SPN2_EXV_POTENTIAL_HPP
#include <mjolnir/forcefield/3SPN2/ThreeSPN2Common.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/math/math.hpp>
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
template<typename realT>
class ThreeSPN2ExcludedVolumePotential
{
  public:
    using real_type = realT;
    using bead_kind = parameter_3SPN2::bead_kind;

    using parameter_type      = bead_kind;
    using pair_parameter_type = real_type; // sigma_ij

    using container_type      = std::vector<parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList;

    static constexpr parameter_type default_parameter() noexcept
    {
        return bead_kind::Unknown;
    }

  public:

    ThreeSPN2ExcludedVolumePotential(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        ignore_group_type ignore_grp)
        : exclusion_list_({{"next_nucl",   1}, {"bond",           1},
                           {"3SPN2_angle", 1}, {"3SPN2_dihedral", 1}},
                          ignore_molecule_type("Nothing"), std::move(ignore_grp))
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
    ~ThreeSPN2ExcludedVolumePotential() = default;

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return (this->sigma_of(this->parameters_[i]) +
                this->sigma_of(this->parameters_[j])) * 0.5;
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

    real_type potential(const real_type r, const pair_parameter_type& sigma) const noexcept
    {
        if(sigma <= r){return 0;}

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return this->epsilon_ * (r12s12 - 2 * r6s6 + 1);
    }
    real_type derivative(const real_type r, const pair_parameter_type& sigma) const noexcept
    {
        if(sigma <= r){return 0;}

        const real_type rinv   = 1 / r;
        const real_type r1s1   = sigma * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 12 * epsilon_ * rinv * (r6s6 - r12s12);
    }

    real_type max_cutoff_length() const noexcept
    {
        return this->cutoff_;
    }

    template<typename traitsT>
    void initialize(const System<traitsT>& sys) noexcept
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
                                 unit::constants<real_type>::J_to_cal);

                // convert from kJ/mol to kcal/mol (/= 4.18)
                this->epsilon_ *= unit::constants<real_type>::J_to_cal;
            }

            const auto& length_unit = physics::constants<real_type>::length_unit();
            MJOLNIR_LOG_INFO("length unit is ", length_unit);
            assert(length_unit == "nm" || length_unit == "angstrom");

            if(length_unit == "nm")
            {
                MJOLNIR_LOG_INFO("length unit (nm) differs from the default, "
                                 "[angstrom]. converting by multiplying ",
                                 unit::constants<real_type>::angstrom_to_nm);

                // convert angstrom -> nm (* 0.1)
                this->cutoff_   *= unit::constants<real_type>::angstrom_to_nm;

                this->sigma_P_  *= unit::constants<real_type>::angstrom_to_nm;
                this->sigma_S_  *= unit::constants<real_type>::angstrom_to_nm;
                this->sigma_A_  *= unit::constants<real_type>::angstrom_to_nm;
                this->sigma_T_  *= unit::constants<real_type>::angstrom_to_nm;
                this->sigma_G_  *= unit::constants<real_type>::angstrom_to_nm;
                this->sigma_C_  *= unit::constants<real_type>::angstrom_to_nm;
            }
            unit_converted_ = true;
        }
        // construct a exclusion list
        this->update(sys);
        return;
    }

    // nothing to do when system parameters change.
    template<typename traitsT>
    void update(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // make exclusion list based on the topology
        exclusion_list_.make(sys);

        // --------------------------------------------------------------------
        // list up beads that are within 3 nucleotides
        const auto& topol = sys.topology();
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
            auto within3 = topol.list_adjacent_within(i, 3, "next_nucl");

            // remove non-base stuff
            within3.erase(std::remove_if(within3.begin(), within3.end(), is_base),
                          within3.end());
            // make it unique
            std::sort(within3.begin(), within3.end());
            within3.erase(std::unique(within3.begin(), within3.end()),
                          within3.end());

            // push the list to the master list
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

    // ------------------------------------------------------------------------
    // to check bases has base-pairing interaction.
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        if(exclusion_list_.is_excluded(i, j))
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

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}

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
    real_type epsilon_   = 1.0; // [kJ/mol]
    real_type sigma_P_   = 4.5; // [angstrom]
    real_type sigma_S_   = 6.2;
    real_type sigma_A_   = 5.4;
    real_type sigma_T_   = 7.1;
    real_type sigma_G_   = 4.9;
    real_type sigma_C_   = 6.4;

    real_type cutoff_    = 7.1; // the maximum sigma_ij

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

} // mjolnir
#endif // MJOLNIR_POTENTIAL_GLOBAL_3SPN2_EXV_POTENTIAL_HPP
