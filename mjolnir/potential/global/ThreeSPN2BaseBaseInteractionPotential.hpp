#ifndef MJOLNIR_POTENTIAL_GLOBAL_3SPN2_BASE_BASE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_3SPN2_BASE_BASE_POTENTIAL_HPP
#include <mjolnir/potential/global/ThreeSPN2Common.hpp>
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

// It calculates a base-pairing energy/force that is a part of 3SPN2 DNA model.
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
// The potential function closely tied to the interaction, so this potential
// will only be used by ThreeSPN2BaseBaseIntearction.
//
// So the interface of this potential is completely different from other global
// potentials. Please be careful.
//
// Note: an identifier starts with a digit is not allowed in C++ standard.
//       see N3337 2.11 for detail. So `3SPN2BaseBaseInteraction` is not a valid name.
template<typename realT>
class ThreeSPN2BaseBaseInteractionPotential
{
  public:
    using real_type = realT;

    struct parameter_type
    {
        base_kind   base;
        std::size_t strand_index; // index in a strand.
        std::size_t S_idx, B3_idx, B5_idx;
    };
    struct pair_parameter_type
    {
        base_pair_kind BP_kind;
        std::size_t    Si;
        std::size_t    Sj;
    };
    using container_type = std::vector<parameter_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;

    static constexpr std::size_t invalid() noexcept
    {
        return std::numeric_limits<std::size_t>::max();
    }

    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{
            base_kind::X, invalid(), invalid(), invalid(), invalid()
        };
    }

  public:

    ThreeSPN2BaseBaseInteractionPotential(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : cutoff_(18.0), cutoff_sq_(18.0*18.0),
          exclusion_list_(std::move(ignore_mol), std::move(ignore_grp),
                          {{"next_nucleotide", 3}})
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
    ~ThreeSPN2BaseBaseInteractionPotential() = default;

    base_pair_kind bp_kind(const base_kind lhs, const base_kind rhs) const noexcept
    {
        switch(lhs)
        {
            case base_kind::A: {assert(rhs == base_kind::T); return base_pair_kind::AT;}
            case base_kind::T: {assert(rhs == base_kind::A); return base_pair_kind::TA;}
            case base_kind::C: {assert(rhs == base_kind::G); return base_pair_kind::CG;}
            case base_kind::G: {assert(rhs == base_kind::C); return base_pair_kind::GC;}
            default:           {assert(false);}
        }
    }
    cross_stack_kind cs5_kind(const base_kind lhs, const base_kind rhs) const noexcept
    {
        assert(lhs != base_kind::X);
        assert(rhs != base_kind::X);
        const auto lhs_u8 = static_cast<std::uint8_t>(lhs);
        const auto rhs_u8 = static_cast<std::uint8_t>(rhs);

        return static_cast<cross_stack_kind>(lhs_u8 << 2 | rhs_u8);
    }
    cross_stack_kind cs3_kind(const base_kind lhs, const base_kind rhs) const noexcept
    {
        assert(lhs != base_kind::X);
        assert(rhs != base_kind::X);
        const auto lhs_u8 = static_cast<std::uint8_t>(lhs);
        const auto rhs_u8 = static_cast<std::uint8_t>(rhs);

        return static_cast<cross_stack_kind>((lhs_u8 << 2 | rhs_u8) + 16u);
    }

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        const auto& para_i = this->parameters_[i];
        const auto& para_j = this->parameters_[j];
        return pair_parameter_type{
            this->bp_kind(para_i.base, para_j.base), para_i.S_idx, para_j.S_idx
        };
    }

    real_type epsilon(const base_pair_kind bp) const noexcept
    {
        switch(bp)
        {
            case base_pair_kind::AT: {return this->epsilon_AT_;}
            case base_pair_kind::TA: {return this->epsilon_AT_;} // symmetric
            case base_pair_kind::GC: {return this->epsilon_GC_;}
            case base_pair_kind::CG: {return this->epsilon_GC_;} // symmetric
            default{assert(false);}
        }
    }
    real_type epsilon(const cross_stack_kind cs) const noexcept
    {
        return epsilon_CS_[static_cast<std::uint8_t>(cs)];
    }

    real_type alpha(const base_pair_kind) const noexcept
    {
        return this->alpha_BP_;
    }
    real_type alpha(const cross_stack_kind) const noexcept
    {
        return this->alpha_CS_;
    }

    real_type K_BP()         const noexcept {return this->K_BP_;}
    real_type K_CS()         const noexcept {return this->K_CS_;}
    real_type pi_over_K_BP() const noexcept {return this->pi_over_K_BP_;}
    real_type pi_over_K_CS() const noexcept {return this->pi_over_K_CS_;}

    real_type r0(const base_pair_kind bp) const noexcept
    {
        switch(bp)
        {
            case base_pair_kind::AT: {return this->r0_AT_;}
            case base_pair_kind::TA: {return this->r0_AT_;} // symmetric
            case base_pair_kind::GC: {return this->r0_GC_;}
            case base_pair_kind::CG: {return this->r0_GC_;} // symmetric
            default{assert(false);}
        }
    }
    real_type r0(const cross_stac_kind cs) const noexcept
    {
        return r0_CS_[static_cast<std::uint8_t>(cs)];
    }
    real_type theta1_0(const base_pair_kind bp) const noexcept
    {
        switch(bp)
        {
            case base_pair_kind::AT: {return this->theta1_0_AT_;}
            case base_pair_kind::TA: {return this->theta1_0_TA_;}
            case base_pair_kind::GC: {return this->theta1_0_GC_;}
            case base_pair_kind::CG: {return this->theta1_0_CG_;}
            default{assert(false);}
        }
    }
    real_type theta2_0(const base_pair_kind bp) const noexcept
    {
        switch(bp)
        {
            case base_pair_kind::AT: {return this->theta2_0_AT_;}
            case base_pair_kind::TA: {return this->theta2_0_TA_;}
            case base_pair_kind::GC: {return this->theta2_0_GC_;}
            case base_pair_kind::CG: {return this->theta2_0_CG_;}
            default{assert(false);}
        }
    }
    real_type theta3_0(const base_pair_kind bp) const noexcept
    {
        switch(bp)
        {
            case base_pair_kind::AT: {return this->theta3_0_AT_;}
            case base_pair_kind::TA: {return this->theta3_0_TA_;}
            case base_pair_kind::GC: {return this->theta3_0_GC_;}
            case base_pair_kind::CG: {return this->theta3_0_CG_;}
            default{assert(false);}
        }
    }
    real_type thetaCS_0(const cross_stac_kind cs) const noexcept
    {
        return theta_CS_0_[static_cast<std::uint8_t>(cs)];
    }
    real_type phi_0(const base_pair_kind bp) const noexcept
    {
        switch(bp)
        {
            case base_pair_kind::AT: {return this->phi_0_AT_;}
            case base_pair_kind::TA: {return this->phi_0_TA_;}
            case base_pair_kind::GC: {return this->phi_0_GC_;}
            case base_pair_kind::CG: {return this->phi_0_CG_;}
            default{assert(false);}
        }
    }

    real_type f(const base_pair_kind,
                const real_type theta, const real_type theta0) const noexcept
    {
        return this->f_impl(K_BP_, pi_over_K_BP_, theta, theta0);
    }
    real_type f(const cross_stack_kind,
                const real_type theta, const real_type theta0) const noexcept
    {
        return this->f_impl(K_CS_, pi_over_K_CS_, theta, theta0);
    }
    real_type df(const base_pair_kind,
                 const real_type theta, const real_type theta0) const noexcept
    {
        return this->df_impl(K_BP_, pi_over_K_BP_, theta, theta0);
    }
    real_type df(const cross_stack_kind,
                 const real_type theta, const real_type theta0) const noexcept
    {
        return this->df_impl(K_CS_, pi_over_K_CS_, theta, theta0);
    }

    real_type cutoff_sq() const noexcept {return this->cutoff_sq_;}
    real_type cutoff()    const noexcept {return this->cutoff_;}

    real_type max_cutoff_length() const noexcept
    {
        return this->cutoff_;
    }

    template<typename traitsT>
    void initialize(const System<traitsT>& sys) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(!unit_converted)
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
                this->epsilon_AT_ *= unit::constants<real_type>::J_to_cal;
                this->epsilon_GC_ *= unit::constants<real_type>::J_to_cal;

                for(auto& epsilon_CS : this->epsilon_CS_)
                {
                    epsilon_CS *= unit::constants<real_type>::J_to_cal;
                }
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
                this->cutoff_sq_ = this->cutoff_ * this->cutoff_;

                this->r0_AT_ *= unit::constants<real_type>::angstrom_to_nm;
                this->r0_GC_ *= unit::constants<real_type>::angstrom_to_nm;

                for(auto& r0_CS : this->r0_CS_)
                {
                    r0_CS *= unit::constants<real_type>::angstrom_to_nm;
                }
            }

            MJOLNIR_LOG_INFO("angle parameters are convered into rad.");
            for(auto& theta_CS_0 : this->theta_CS_0_)
            {
                theta_CS_0 *= (math::constants<real_type>::pi / 180.0);
            }
            unit_converted = true;
        }
        // construct a exclusion list
        return;
    }

    // nothing to do when system parameters change.
    template<typename traitsT>
    void update(const System<traitsT>&) noexcept {return;}

    // ------------------------------------------------------------------------
    // to check bases has base-pairing interaction.
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        switch(this->parameters_[i].base)
        {
            case base_kind::A: {return this->parameters_[j].base == base_kind::T;}
            case base_kind::T: {return this->parameters_[j].base == base_kind::A;}
            case base_kind::C: {return this->parameters_[j].base == base_kind::G;}
            case base_kind::G: {return this->parameters_[j].base == base_kind::C;}
            default:           {assert(false);}
        }
    }

    exclusion_list_type const& exclusion_list() const noexcept {return exclusion_list_;}

    // ------------------------------------------------------------------------
    // used by Observer.
    static const char* name() noexcept {return "ThreeSPN2BaseBaseIntearction";}

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    std::vector<std::size_t> const& participants() const noexcept {return participants_;}

  private:

    real_type f_impl(const real_type K,     const real_type pi_over_K,
                     const real_type theta, const real_type theta0) const noexcept
    {
        const auto dtheta     = theta - theta0;
        const auto abs_dtheta = std::abs(dtheta);
        if(abs_dtheta < pi_over_K * 0.5)
        {
            return 1.0;
        }
        else if(abs_dtheta < pi_over_K)
        {
            const auto cos_Kdtheta = std::cos(K * dtheta);
            return 1.0 - cos_Kdtheta * cos_Kdtheta;
        }
        else
        {
            return 0.0;
        }
    }
    real_type df_impl(const real_type K,     const real_type pi_over_K,
                      const real_type theta, const real_type theta0) const noexcept
    {
        const auto dtheta     = theta - theta0;
        const auto abs_dtheta = std::abs(dtheta);

        if(abs_dtheta < pi_over_K * 0.5)
        {
            return 0.0;
        }
        else if(abs_dtheta < pi_over_K)
        {
            return K * std::sin(2 * K * dtheta);
        }
        else
        {
            return 0.0;
        }
    }

  private:

    // the default values are in Angstrom and kJ/mol, as in the original paper.
    // in `initialize()`, the units will be converted.
    bool unit_converted     = false;

    real_type cutoff_       = 18.0;        // [angstrom]
    real_type cutoff_sq_    = 18.0 * 18.0; // [angstrom^2]

    real_type alpha_BP_     =  2.0;
    real_type alpha_CS_     =  4.0;

    real_type K_BP_         = 12.0;
    real_type K_CS_         =  8.0;
    real_type pi_over_K_BP_ = math::constants<real_type>::pi / 12.0;
    real_type pi_over_K_BP_ = math::constants<real_type>::pi /  8.0;

    real_type epsilon_AT_   = 16.73; // [kJ/mol]
    real_type epsilon_GC_   = 21.18; // [kJ/mol]

    real_type r0_AT_        = 5.941; // [angstrom]
    real_type r0_GC_        = 5.518; // [angstrom]

    real_type theta1_0_AT_  = 156.54 / 180.0 * math::constants<real_type>::pi;
    real_type theta1_0_TA_  = 135.78 / 180.0 * math::constants<real_type>::pi;
    real_type theta1_0_GC_  = 154.62 / 180.0 * math::constants<real_type>::pi;
    real_type theta1_0_CG_  = 152.74 / 180.0 * math::constants<real_type>::pi;

    real_type theta2_0_AT_  = 135.78 / 180.0 * math::constants<real_type>::pi;
    real_type theta2_0_TA_  = 156.54 / 180.0 * math::constants<real_type>::pi;
    real_type theta2_0_GC_  = 152.74 / 180.0 * math::constants<real_type>::pi;
    real_type theta2_0_CG_  = 154.62 / 180.0 * math::constants<real_type>::pi;

    real_type theta3_0_AT_  = 116.09 / 180.0 * math::constants<real_type>::pi;
    real_type theta3_0_TA_  = 116.09 / 180.0 * math::constants<real_type>::pi;
    real_type theta3_0_GC_  = 131.78 / 180.0 * math::constants<real_type>::pi;
    real_type theta3_0_CG_  = 131.78 / 180.0 * math::constants<real_type>::pi;

    real_type phi_0_AT_     = -38.35 / 180.0 * math::constants<real_type>::pi;
    real_type phi_0_TA_     = -38.35 / 180.0 * math::constants<real_type>::pi;
    real_type phi_0_GC_     = -45.81 / 180.0 * math::constants<real_type>::pi;
    real_type phi_0_CG_     = -45.81 / 180.0 * math::constants<real_type>::pi;

    std::array<real_type, 32> epsilon_CS_ = {{ // [kJ/mol]
        /* AA5 */ 2.186, /* AT5 */ 2.774, /* AG5 */ 2.833, /* AC5 */ 1.951,
        /* TA5 */ 2.774, /* TT5 */ 2.186, /* TG5 */ 2.539, /* TC5 */ 2.980,
        /* GA5 */ 2.833, /* GT5 */ 2.539, /* GG5 */ 3.774, /* GC5 */ 1.129,
        /* CA5 */ 1.951, /* CT5 */ 2.980, /* CG5 */ 1.129, /* CC5 */ 4.802,

        /* AA3 */ 2.186, /* AT3 */ 2.774, /* AG3 */ 2.980, /* AC3 */ 2.539,
        /* TA3 */ 2.774, /* TT3 */ 2.186, /* TG3 */ 1.951, /* TC3 */ 2.833,
        /* GA3 */ 2.980, /* GT3 */ 1.951, /* GG3 */ 4.802, /* GC3 */ 1.129,
        /* CA3 */ 2.539, /* CT3 */ 2.833, /* CG3 */ 1.129, /* CC3 */ 3.774
    }};
    std::array<real_type, 32> r0_CS_      = {{ // [angstrom]
        /* AA5 */ 6.208, /* AT5 */ 6.876, /* AG5 */ 6.072, /* AC5 */ 6.941,
        /* TA5 */ 6.876, /* TT5 */ 7.480, /* TG5 */ 6.771, /* TC5 */ 7.640,
        /* GA5 */ 6.072, /* GT5 */ 6.771, /* GG5 */ 5.921, /* GC5 */ 6.792,
        /* CA5 */ 6.941, /* CT5 */ 7.640, /* CG5 */ 6.792, /* CC5 */ 7.698,

        /* AA3 */ 5.435, /* AT3 */ 6.295, /* AG3 */ 5.183, /* AC3 */ 5.965,
        /* TA3 */ 6.295, /* TT3 */ 7.195, /* TG3 */ 6.028, /* TC3 */ 6.868,
        /* GA3 */ 5.183, /* GT3 */ 6.028, /* GG3 */ 4.934, /* GC3 */ 5.684,
        /* CA3 */ 5.965, /* CT3 */ 6.868, /* CG3 */ 5.684, /* CC3 */ 6.453
    }};
    // XXX in [degree]. would be converted to [radian] in `initialize()`
    std::array<real_type, 32> theta_CS_0_ = {{
        /* AA5 */ 154.38, /* AT5 */ 159.10, /* AG5 */ 152.46, /* AC5 */ 157.58,
        /* TA5 */ 147.10, /* TT5 */ 153.79, /* TG5 */ 144.44, /* TC5 */ 148.59,
        /* GA5 */ 154.69, /* GT5 */ 157.83, /* GG5 */ 153.43, /* GC5 */ 158.60,
        /* CA5 */ 160.37, /* CT5 */ 164.45, /* CG5 */ 158.62, /* CC5 */ 162.73,

        /* AA3 */ 116.88, /* AT3 */ 121.74, /* AG3 */ 114.23, /* AC3 */ 114.58,
        /* TA3 */ 109.42, /* TT3 */ 112.95, /* TG3 */ 107.32, /* TC3 */ 106.41,
        /* GA3 */ 119.34, /* GT3 */ 124.72, /* GG3 */ 116.51, /* GC3 */ 117.49,
        /* CA3 */ 122.10, /* CT3 */ 125.80, /* CG3 */ 120.00, /* CC3 */ 119.67
    }};

    container_type           parameters_;   // indices
    std::vector<std::size_t> participants_; // base beads

    exclusion_list_type exclusion_list_;
};

} // mjolnir
#endif // MJOLNIR_POTENTIAL_GLOBAL_3SPN2_BASE_BASE_POTENTIAL_HPP
