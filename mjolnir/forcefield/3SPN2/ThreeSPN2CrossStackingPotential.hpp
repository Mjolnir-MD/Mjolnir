#ifndef MJOLNIR_FORCEFIELD_3SPN2_CROSS_STACKING_POTENTIAL_HPP
#define MJOLNIR_FORCEFIELD_3SPN2_CROSS_STACKING_POTENTIAL_HPP
#include <mjolnir/forcefield/global/ParameterList.hpp>
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

// It calculates cross-stacking energy/force that is a part of 3SPN2 DNA model.
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
// And also 3SPN2.C model.
// - G. S. Freeman, D. M. Hinckley, J. P. Lequieu, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2014)
//
// The potential function closely tied to the interaction, so this potential
// will only be used by ThreeSPN2CrossStackingInteraction.
//
// So the interface of this potential is completely different from other global
// potentials. Please be careful.
//
// Note: an identifier starts with a digit is not allowed in C++ standard.
//       see N3337 2.11 for detail. So `3SPN2BaseBaseInteraction` is not a valid name.

//       Si   Bi   Bj   Sj
//  5'    o -- o===o -- o     3'
//  ^    /      \        \    |
//  | P o        \        o P |
//  |    \        \      /    v
//  3'    o -- o===o -- o     5'
//       Bi_next   Bj_next
//
template<typename realT>
class ThreeSPN2CrossStackingPotential
{
  public:
    using real_type        = realT;
    using base_kind        = parameter_3SPN2::base_kind;
    using base_pair_kind   = parameter_3SPN2::base_pair_kind;
    using cross_stack_kind = parameter_3SPN2::cross_stack_kind;

    struct parameter_type // per pair (Bi, Bj)
    {
        base_pair_kind bp_kind; // Base pair kind of the pairing baeds
    };

    static constexpr std::size_t invalid() noexcept
    {
        return std::numeric_limits<std::size_t>::max();
    }

  public:

    ThreeSPN2CrossStackingPotential() noexcept {}

    real_type f(const real_type K,     const real_type pi_over_K,
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

    real_type df(const real_type K,     const real_type pi_over_K,
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

    real_type U_attr(const real_type epsilon, const real_type alpha,
                     const real_type r,       const real_type r0) const noexcept
    {
        // --------------------------------------------------------
        // U_m^attr =
        //   -e                             ... (dBij <= dBij0)
        //   -e + e * (1 - exp(-a(r-r0)))^2 ... (otherwise)
        //
        if(r <= r0)
        {
            return -epsilon;
        }
        else
        {
            const auto term = real_type(1) - std::exp(-alpha * (r - r0));
            return epsilon * (term * term - real_type(1));
        }
    }

    real_type dU_attr(const real_type epsilon, const real_type alpha,
                      const real_type r,       const real_type r0) const noexcept
    {
        // --------------------------------------------------------
        // dU_m^attr / dr =
        //   0                                 ... (dBij <= dBij0)
        //   2ae(1-exp(-a(r-r0)))exp(-a(r-r0)) ... (otherwise)
        //
        if(r <= r0) {return real_type(0);}
        const auto term = std::exp(-alpha * (r - r0));
        return 2 * alpha * epsilon * term * (real_type(1) - term);
    }

    std::pair<real_type, real_type> U_dU_attr(
            const real_type epsilon, const real_type alpha,
            const real_type r,       const real_type r0) const noexcept
    {
        if(r <= r0)
        {
            return std::make_pair(-epsilon, 0.0);
        }

        const auto term1 = std::exp(-alpha * (r - r0));
        const auto term2 = real_type(1) - term1;
        return std::make_pair(epsilon * (term2 * term2 - real_type(1)),
                              2 * alpha * epsilon * term1 * term2);
    }

    real_type U_rep(const real_type epsilon, const real_type alpha,
                    const real_type r,       const real_type r0) const noexcept
    {
        if(r0 < r) {return real_type(0);}
        const auto term = real_type(1) - std::exp(-alpha * (r - r0));
        return epsilon * term * term;
    }
    real_type dU_rep(const real_type epsilon, const real_type alpha,
                     const real_type r,       const real_type r0) const noexcept
    {
        if(r0 < r) {return real_type(0);}
        const auto term = std::exp(-alpha * (r - r0));
        return 2 * alpha * epsilon * term * (real_type(1) - term);
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    static const char* name() noexcept {return "3SPN2CrossStacking";}
};

template<typename traitsT>
class ThreeSPN2CrossStackingParameterList
    : public ParameterListBase<traitsT,
        ThreeSPN2CrossStackingPotential<typename traitsT::real_type>>
{
  public:
    using traits_type         = traitsT;
    using real_type           = typename traits_type::real_type;
    using potential_type      = ThreeSPN2CrossStackingPotential<real_type>;
    using base_type           = ParameterListBase<traits_type, potential_type>;
    using pair_parameter_type = typename base_type::pair_parameter_type;
    // pair_parameter_type = potential_type::parameter_type;

    // topology stuff
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using molecule_id_type     = typename base_type::molecule_id_type;
    using group_id_type        = typename base_type::group_id_type;
    using connection_kind_type = typename base_type::connection_kind_type;
    using ignore_molecule_type = typename base_type::ignore_molecule_type;
    using ignore_group_type    = typename base_type::ignore_group_type;
    using exclusion_list_type  = typename base_type::exclusion_list_type;

    using base_kind        = parameter_3SPN2::base_kind;
    using base_pair_kind   = parameter_3SPN2::base_pair_kind;
    using cross_stack_kind = parameter_3SPN2::cross_stack_kind;

    struct parameter_type // per particle
    {
        base_kind   base;
        std::size_t S, B3, B5;
    };
    using container_type = std::vector<parameter_type>;

  public:

    template<typename ParameterSet>
    ThreeSPN2CrossStackingParameterList(ParameterSet para_set,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>&         exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : cutoff_      (para_set.cutoff),
          cutoff_sq_   (para_set.cutoff_sq),
          alpha_CS_    (para_set.alpha_CS),
          K_CS_        (para_set.K_CS),
          K_BP_        (para_set.K_BP),
          pi_over_K_CS_(para_set.pi_over_K_CS),
          pi_over_K_BP_(para_set.pi_over_K_BP),
          // --------------------------------
          theta3_0_AT_(para_set.theta3_0_AT),
          theta3_0_TA_(para_set.theta3_0_TA),
          theta3_0_GC_(para_set.theta3_0_GC),
          theta3_0_CG_(para_set.theta3_0_CG),
          // --------------------------------
          epsilon_CS_(para_set.epsilon_CS),
          r0_CS_     (para_set.r0_CS),
          theta_CS_0_(para_set.theta_CS_0),
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

        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, parameter_type{base_kind::X,
                        potential_type::invalid(), potential_type::invalid(),
                        potential_type::invalid()
                    });
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~ThreeSPN2CrossStackingParameterList() override {}
    ThreeSPN2CrossStackingParameterList(const ThreeSPN2CrossStackingParameterList&) = default;
    ThreeSPN2CrossStackingParameterList(ThreeSPN2CrossStackingParameterList&&)      = default;
    ThreeSPN2CrossStackingParameterList& operator=(const ThreeSPN2CrossStackingParameterList&) = default;
    ThreeSPN2CrossStackingParameterList& operator=(ThreeSPN2CrossStackingParameterList&&)      = default;

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

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        const auto& para_i_base = this->parameters_[i].base;
        const auto& para_j_base = this->parameters_[j].base;
        return pair_parameter_type{this->bp_kind(para_i_base, para_j_base)};
    }

    real_type epsilon(const cross_stack_kind cs) const noexcept
    {
        return epsilon_CS_[static_cast<std::uint8_t>(cs)];
    }

    real_type alpha(const cross_stack_kind) const noexcept
    {
        return this->alpha_CS_;
    }
    real_type K_CS()         const noexcept {return this->K_CS_;}
    real_type K_BP()         const noexcept {return this->K_BP_;}
    real_type pi_over_K_CS() const noexcept {return this->pi_over_K_CS_;}
    real_type pi_over_K_BP() const noexcept {return this->pi_over_K_BP_;}

    real_type r0(const cross_stack_kind cs) const noexcept
    {
        return r0_CS_[static_cast<std::uint8_t>(cs)];
    }
    real_type theta3_0(const base_pair_kind bp) const noexcept
    {
        switch(bp)
        {
            case base_pair_kind::AT: {return this->theta3_0_AT_;}
            case base_pair_kind::TA: {return this->theta3_0_TA_;}
            case base_pair_kind::GC: {return this->theta3_0_GC_;}
            case base_pair_kind::CG: {return this->theta3_0_CG_;}
            default: {assert(false);}
        }
    }
    real_type thetaCS_0(const cross_stack_kind cs) const noexcept
    {
        return theta_CS_0_[static_cast<std::uint8_t>(cs)];
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
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
                                 unit::constants<real_type>::J_to_cal());

                // convert from kJ/mol to kcal/mol (/= 4.18)
                for(auto& epsilon_CS : this->epsilon_CS_)
                {
                    epsilon_CS *= unit::constants<real_type>::J_to_cal();
                }
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
                this->cutoff_sq_ = this->cutoff_ * this->cutoff_;

                for(auto& r0_CS : this->r0_CS_)
                {
                    r0_CS *= unit::constants<real_type>::angstrom_to_nm();
                }
            }

            MJOLNIR_LOG_INFO("angle parameters are convered into rad.");
            for(auto& theta_CS_0 : this->theta_CS_0_)
            {
                theta_CS_0 *= (math::constants<real_type>::pi() / 180.0);
            }
            unit_converted = true;
        }
        // construct a exclusion list
        this->update(sys, topol, pot);
        return;
    }

    // nothing to do when system parameters change.
    void update(const system_type& sys, const topology_type& topol,
                const potential_type&) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        exclusion_list_.make(sys, topol);
        return;
    }

    // -----------------------------------------------------------------------
    // for spatial partitions

    real_type max_cutoff_length() const noexcept override
    {
        return this->cutoff_;
    }

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

    // to check bases has base-pairing interaction.
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        if(j <= i || exclusion_list_.is_excluded(i, j))
        {
            return false;
        }
        switch(this->parameters_[i].base)
        {
            case base_kind::A: {return this->parameters_[j].base == base_kind::T;}
            case base_kind::T: {return this->parameters_[j].base == base_kind::A;}
            case base_kind::C: {return this->parameters_[j].base == base_kind::G;}
            case base_kind::G: {return this->parameters_[j].base == base_kind::C;}
            default:           {assert(false);}
        }
    }

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    base_type* clone() const override
    {
        return new ThreeSPN2CrossStackingParameterList(*this);
    }

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return parameters_;}
    std::vector<parameter_type> const& parameters() const noexcept {return parameters_;}

    parameter_type const& at(std::size_t i) const noexcept {return parameters_.at(i);}
    parameter_type&       at(std::size_t i)       noexcept {return parameters_.at(i);}

    real_type cutoff_sq() const noexcept {return this->cutoff_sq_;}
    real_type cutoff()    const noexcept {return this->cutoff_;}


  private:

    // the default values are in Angstrom and kJ/mol, as in the original paper.
    // in `initialize()`, the units will be converted.
    bool unit_converted     = false;

    // -----------------------------------------------------------------------
    // forcefield parameters

    real_type cutoff_;    // initially, [angstrom]
    real_type cutoff_sq_; //            [angstrom^2]

    real_type alpha_CS_;
    real_type K_CS_;
    real_type K_BP_;
    real_type pi_over_K_CS_;
    real_type pi_over_K_BP_;

    real_type theta3_0_AT_;
    real_type theta3_0_TA_;
    real_type theta3_0_GC_;
    real_type theta3_0_CG_;

    std::array<real_type, 32> epsilon_CS_; // [kJ/mol]
    std::array<real_type, 32> r0_CS_     ; // [angstrom]
    std::array<real_type, 32> theta_CS_0_; // [degree] XXX should be converted in initialize()

    // -----------------------------------------------------------------------
    // runtime parameters

    container_type           parameters_;   // indices
    std::vector<std::size_t> participants_; // base beads

    exclusion_list_type exclusion_list_;
};

// parameter set for 3SPN2.
template<typename realT>
struct ThreeSPN2CrossStackingGlobalPotentialParameter
{
    using real_type = realT;

    real_type cutoff       = 18.0;        // [angstrom]
    real_type cutoff_sq    = 18.0 * 18.0; // [angstrom^2]

    real_type alpha_CS     =  4.0;
    real_type K_CS         =  8.0;
    real_type K_BP         = 12.0;
    real_type pi_over_K_CS = math::constants<real_type>::pi() /  8.0;
    real_type pi_over_K_BP = math::constants<real_type>::pi() / 12.0;

    real_type theta3_0_AT  = 116.09 / 180.0 * math::constants<real_type>::pi();
    real_type theta3_0_TA  = 116.09 / 180.0 * math::constants<real_type>::pi();
    real_type theta3_0_GC  = 124.94 / 180.0 * math::constants<real_type>::pi();
    real_type theta3_0_CG  = 124.94 / 180.0 * math::constants<real_type>::pi();

    std::array<real_type, 32> epsilon_CS = {{ // [kJ/mol]
        /* AA5 */ 2.186, /* AT5 */ 2.774, /* AG5 */ 2.833, /* AC5 */ 1.951,
        /* TA5 */ 2.774, /* TT5 */ 2.186, /* TG5 */ 2.539, /* TC5 */ 2.980,
        /* GA5 */ 2.833, /* GT5 */ 2.539, /* GG5 */ 3.774, /* GC5 */ 1.129,
        /* CA5 */ 1.951, /* CT5 */ 2.980, /* CG5 */ 1.129, /* CC5 */ 4.802,

        /* AA3 */ 2.186, /* AT3 */ 2.774, /* AG3 */ 2.980, /* AC3 */ 2.539,
        /* TA3 */ 2.774, /* TT3 */ 2.186, /* TG3 */ 1.951, /* TC3 */ 2.833,
        /* GA3 */ 2.980, /* GT3 */ 1.951, /* GG3 */ 4.802, /* GC3 */ 1.129,
        /* CA3 */ 2.539, /* CT3 */ 2.833, /* CG3 */ 1.129, /* CC3 */ 3.774
    }};
    std::array<real_type, 32> r0_CS      = {{ // [angstrom]
        /* AA5 */ 6.208, /* AT5 */ 6.876, /* AG5 */ 6.072, /* AC5 */ 6.811,
        /* TA5 */ 6.876, /* TT5 */ 7.480, /* TG5 */ 6.771, /* TC5 */ 7.453,
        /* GA5 */ 6.072, /* GT5 */ 6.771, /* GG5 */ 5.921, /* GC5 */ 6.688,
        /* CA5 */ 6.811, /* CT5 */ 7.453, /* CG5 */ 6.688, /* CC5 */ 7.409,

        /* AA3 */ 5.435, /* AT3 */ 6.295, /* AG3 */ 5.183, /* AC3 */ 6.082,
        /* TA3 */ 6.295, /* TT3 */ 7.195, /* TG3 */ 6.028, /* TC3 */ 6.981,
        /* GA3 */ 5.183, /* GT3 */ 6.028, /* GG3 */ 4.934, /* GC3 */ 5.811,
        /* CA3 */ 6.082, /* CT3 */ 6.981, /* CG3 */ 5.811, /* CC3 */ 6.757
    }};
    // XXX in [degree]. would be converted to [radian] in `initialize()`
    std::array<real_type, 32> theta_CS_0 = {{
        /* AA5 */ 154.38, /* AT5 */ 159.10, /* AG5 */ 152.46, /* AC5 */ 158.38,
        /* TA5 */ 147.10, /* TT5 */ 153.79, /* TG5 */ 144.44, /* TC5 */ 151.48,
        /* GA5 */ 154.69, /* GT5 */ 157.83, /* GG5 */ 153.43, /* GC5 */ 158.04,
        /* CA5 */ 152.99, /* CT5 */ 159.08, /* CG5 */ 150.53, /* CC5 */ 157.17,

        /* AA3 */ 116.88, /* AT3 */ 121.74, /* AG3 */ 114.23, /* AC3 */ 119.06,
        /* TA3 */ 109.42, /* TT3 */ 112.95, /* TG3 */ 107.32, /* TC3 */ 110.56,
        /* GA3 */ 119.34, /* GT3 */ 124.72, /* GG3 */ 116.51, /* GC3 */ 121.98,
        /* CA3 */ 114.60, /* CT3 */ 118.26, /* CG3 */ 112.45, /* CC3 */ 115.88
    }};
};

// parameter for 3SPN2.C
template<typename realT>
struct ThreeSPN2CCrossStackingGlobalPotentialParameter
{
    using real_type = realT;

    real_type cutoff       = 18.0;        // [angstrom]
    real_type cutoff_sq    = 18.0 * 18.0; // [angstrom^2]

    real_type alpha_CS     =  4.0;
    real_type K_CS         =  8.0;
    real_type K_BP         = 12.0;
    real_type pi_over_K_CS = math::constants<real_type>::pi() /  8.0;
    real_type pi_over_K_BP = math::constants<real_type>::pi() / 12.0;

    real_type theta3_0_AT  = 110.92 / 180.0 * math::constants<real_type>::pi();
    real_type theta3_0_TA  = 110.92 / 180.0 * math::constants<real_type>::pi();
    real_type theta3_0_GC  = 120.45 / 180.0 * math::constants<real_type>::pi();
    real_type theta3_0_CG  = 120.45 / 180.0 * math::constants<real_type>::pi();

    std::array<real_type, 32> epsilon_CS = {{ // [kJ/mol]
        /* AA5 */ 1.882, /* AT5 */ 2.388, /* AG5 */ 2.439, /* AC5 */ 1.680,
        /* TA5 */ 2.388, /* TT5 */ 1.882, /* TG5 */ 2.187, /* TC5 */ 2.566,
        /* GA5 */ 2.439, /* GT5 */ 2.187, /* GG5 */ 3.250, /* GC5 */ 0.972,
        /* CA5 */ 1.680, /* CT5 */ 2.566, /* CG5 */ 0.972, /* CC5 */ 4.135,

        /* AA3 */ 1.882, /* AT3 */ 2.388, /* AG3 */ 2.566, /* AC3 */ 2.187,
        /* TA3 */ 2.388, /* TT3 */ 1.882, /* TG3 */ 1.680, /* TC3 */ 2.439,
        /* GA3 */ 2.566, /* GT3 */ 1.680, /* GG3 */ 4.135, /* GC3 */ 0.972,
        /* CA3 */ 2.187, /* CT3 */ 2.439, /* CG3 */ 0.972, /* CC3 */ 3.250
    }};

    std::array<real_type, 32> r0_CS      = {{ // [angstrom]
        /* AA5 */ 6.420, /* AT5 */ 6.770, /* AG5 */ 6.270, /* AC5 */ 6.840,
        /* TA5 */ 6.770, /* TT5 */ 7.210, /* TG5 */ 6.530, /* TC5 */ 7.080,
        /* GA5 */ 6.270, /* GT5 */ 6.530, /* GG5 */ 5.740, /* GC5 */ 6.860,
        /* CA5 */ 6.840, /* CT5 */ 7.080, /* CG5 */ 6.860, /* CC5 */ 6.790,

        /* AA3 */ 5.580, /* AT3 */ 6.140, /* AG3 */ 5.630, /* AC3 */ 6.180,
        /* TA3 */ 6.140, /* TT3 */ 6.800, /* TG3 */ 6.070, /* TC3 */ 6.640,
        /* GA3 */ 5.630, /* GT3 */ 6.070, /* GG3 */ 5.870, /* GC3 */ 5.660,
        /* CA3 */ 6.180, /* CT3 */ 6.640, /* CG3 */ 5.660, /* CC3 */ 6.800
    }};

    // XXX in [degree]. would be converted to [radian] in `initialize()`
    std::array<real_type, 32> theta_CS_0 = {{
        /* AA5 */ 154.04, /* AT5 */ 158.77, /* AG5 */ 153.88, /* AC5 */ 157.69,
        /* TA5 */ 148.62, /* TT5 */ 155.05, /* TG5 */ 147.54, /* TC5 */ 153.61,
        /* GA5 */ 153.91, /* GT5 */ 155.72, /* GG5 */ 151.84, /* GC5 */ 157.80,
        /* CA5 */ 152.04, /* CT5 */ 157.72, /* CG5 */ 151.65, /* CC5 */ 154.49,

        /* AA3 */ 116.34, /* AT3 */ 119.61, /* AG3 */ 115.19, /* AC3 */ 120.92,
        /* TA3 */ 107.40, /* TT3 */ 110.76, /* TG3 */ 106.33, /* TC3 */ 111.57,
        /* GA3 */ 121.61, /* GT3 */ 124.92, /* GG3 */ 120.52, /* GC3 */ 124.88,
        /* CA3 */ 112.45, /* CT3 */ 115.43, /* CG3 */ 110.51, /* CC3 */ 115.80
    }};
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize major use-cases
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{
extern template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ThreeSPN2CrossStackingParameterList<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif// MJOLNIR_SEPARATE_BUILD

#endif // MJOLNIR_POTENTIAL_GLOBAL_3SPN2_BASE_BASE_POTENTIAL_HPP
