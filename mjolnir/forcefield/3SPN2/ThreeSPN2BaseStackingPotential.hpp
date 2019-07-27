#ifndef MJOLNIR_POTENTIAL_LOCAL_3SPN2_BASE_STACKING_INTEARACTION_HPP
#define MJOLNIR_POTENTIAL_LOCAL_3SPN2_BASE_STACKING_INTEARACTION_HPP
#include <mjolnir/forcefield/3SPN2/ThreeSPN2Common.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/math/math.hpp>

namespace mjolnir
{

// It calculates a stacking energy/force that is a part of 3SPN2 DNA model.
// - D. M. Hinckley, G. S. Freeman, J. K. Whitmer, and J. J. de Pablo
//   J. Chem. Phys. (2013)
//
// The potential function closely tied to the interaction, so this interaction
// class is implemented with its own, uninterchangeable potential class.
// It does not take any potential class as a template parameter because
// exchanging potential function of this interaction does not make any sense.
//
// Note: an identifier starts with a digit is not allowed in C++ standard.
//       see N3337 2.11 for detail. So `3SPN2BaseStacking` is not a valid name.
template<typename realT>
class ThreeSPN2BaseStackingPotential
{
  public:
    using real_type       = realT;
    using base_kind       = parameter_3SPN2::base_kind;
    using base_stack_kind = parameter_3SPN2::base_stack_kind;
    using parameter_type  = base_stack_kind;

  public:

    ThreeSPN2BaseStackingPotential() = default;
    ~ThreeSPN2BaseStackingPotential() = default;

    ThreeSPN2BaseStackingPotential(const ThreeSPN2BaseStackingPotential&) = default;
    ThreeSPN2BaseStackingPotential(ThreeSPN2BaseStackingPotential&&)      = default;
    ThreeSPN2BaseStackingPotential& operator=(const ThreeSPN2BaseStackingPotential&) = default;
    ThreeSPN2BaseStackingPotential& operator=(ThreeSPN2BaseStackingPotential&&)      = default;

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
                std::uint8_t idx = 0;
                for(auto& epsilon : this->epsilon_BS_)
                {
                    const auto bs = static_cast<base_stack_kind>(idx);
                    MJOLNIR_LOG_INFO_NO_LF("epsilon:", bs, " = ", epsilon,
                                           " [kJ/mol] -> ");
                    epsilon *= unit::constants<real_type>::J_to_cal;
                    MJOLNIR_LOG_INFO(epsilon, "[kcal/mol]");
                    ++idx;
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
                std::uint8_t idx = 0;
                for(auto& r0_BS : this->r0_BS_)
                {
                    const auto bs = static_cast<base_stack_kind>(idx);
                    MJOLNIR_LOG_INFO_NO_LF("r0:", bs, " = ", r0_BS,
                                           " [angstrom] -> ");
                    r0_BS *= unit::constants<real_type>::angstrom_to_nm;
                    MJOLNIR_LOG_INFO(r0_BS, "[nm]");
                    ++idx;
                }
            }

            MJOLNIR_LOG_INFO("angle parameters are convered into rad.");
            std::uint8_t idx = 0;
            for(auto& theta0_BS : this->theta0_BS_)
            {
                const auto bs = static_cast<base_stack_kind>(idx);
                MJOLNIR_LOG_INFO_NO_LF("theta0:", bs, " = ", theta0_BS, " [deg] -> ");
                theta0_BS *= (math::constants<real_type>::pi / 180.0);
                MJOLNIR_LOG_INFO(theta0_BS, "[rad]");
                ++idx;
            }
            unit_converted_ = true;
        }
        this->update(sys);
        return;
    }

    template<typename traitsT>
    void update(const System<traitsT>&) noexcept {}

    base_stack_kind bs_kind(const base_kind lhs, const base_kind rhs) const noexcept
    {
        assert(lhs != base_kind::X);
        assert(rhs != base_kind::X);
        const auto lhs_u8 = static_cast<std::uint8_t>(lhs);
        const auto rhs_u8 = static_cast<std::uint8_t>(rhs);

        return static_cast<base_stack_kind>(lhs_u8 << 2 | rhs_u8);
    }

    real_type r0(const base_stack_kind bs) const noexcept
    {
        return r0_BS_[static_cast<std::uint8_t>(bs)];
    }
    real_type theta_0(const base_stack_kind bs) const noexcept
    {
        return theta0_BS_[static_cast<std::uint8_t>(bs)];
    }

    real_type epsilon(const base_stack_kind bs) const noexcept
    {
        return epsilon_BS_[static_cast<std::uint8_t>(bs)];
    }
    real_type alpha(const base_stack_kind) const noexcept
    {
        return this->alpha_BS_;
    }

    real_type K_BS()         const noexcept {return this->K_BS_;}
    real_type pi_over_K_BS() const noexcept {return this->pi_over_K_BS_;}

    real_type f(const real_type theta, const real_type theta0) const noexcept
    {
        const auto dtheta     = theta - theta0;
        const auto abs_dtheta = std::abs(dtheta);
        if(abs_dtheta < this->pi_over_K_BS_ * 0.5)
        {
            return 1.0;
        }
        else if(abs_dtheta < this->pi_over_K_BS_)
        {
            const auto cos_Kdtheta = std::cos(this->K_BS_ * dtheta);
            return 1.0 - cos_Kdtheta * cos_Kdtheta;
        }
        else
        {
            return 0.0;
        }
    }
    real_type df(const real_type theta, const real_type theta0) const noexcept
    {
        const auto dtheta     = theta - theta0;
        const auto abs_dtheta = std::abs(dtheta);

        if(abs_dtheta < this->pi_over_K_BS_ * 0.5)
        {
            return 0.0;
        }
        else if(abs_dtheta < this->pi_over_K_BS_)
        {
            return K_BS_ * std::sin(2 * K_BS_ * dtheta);
        }
        else
        {
            return 0.0;
        }
    }

    real_type U_rep(const base_stack_kind bs, const real_type r) const noexcept
    {
        // --------------------------------------------------------------------
        // U_rep = e_ij (1 - exp(-a_ij (r - r0)))^2 ... r < r0
        //       = 0                                ... r0 <= r

        const auto r0  = this->r0(bs);
        if(r0 <= r) {return 0.0;}

        const auto eps  = this->epsilon(bs);
        const auto alp  = this->alpha  (bs);
        const auto term = real_type(1.0) - std::exp(-alp * (r - r0));
        return eps * term * term;
    }
    real_type dU_rep(const base_stack_kind bs, const real_type r) const noexcept
    {
        // --------------------------------------------------------------------
        // dU_rep = 2 a e exp(-a(r-r0)) (1-exp(-a(r-r0))) ... r  <  r0
        //        = 0                                     ... r0 <= r

        const auto r0 = this->r0(bs);
        if(r0 <= r) {return 0.0;}

        const auto eps  = this->epsilon(bs);
        const auto alp  = this->alpha  (bs);
        const auto term = std::exp(-alp * (r - r0));

        return 2 * alp * eps * term * (real_type(1) - term);
    }

    real_type U_attr(const base_stack_kind bs, const real_type r) const noexcept
    {
        // --------------------------------------------------------------------
        // U_attr = -e                             ... (r <= r0)
        //          -e + e * (1 - exp(-a(r-r0)))^2 ... (otherwise)
        //
        const auto r0  = this->r0(bs);
        const auto eps = this->epsilon(bs);
        if(r <= r0) {return -eps;}

        const auto alp  = this->alpha(bs);
        const auto term = std::exp(-alp * (r - r0));
        return eps * (term * term - 2 * term);
    }
    real_type dU_attr(const base_stack_kind bs, const real_type r) const noexcept
    {
        // --------------------------------------------------------------------
        // dU_attr =  0                                 ... (r <= r0)
        //            2ae(1-exp(-a(r-r0)))exp(-a(r-r0)) ... (otherwise)

        const auto r0 = this->r0(bs);
        if(r <= r0) {return 0.0;}

        const auto eps  = this->epsilon(bs);
        const auto alp  = this->alpha  (bs);
        const auto term = std::exp(-alp * (r - r0));

        return 2 * alp * eps * term * (real_type(1) - term);
    }
    std::pair<real_type, real_type>
    U_dU_attr(const base_stack_kind bs, const real_type r) const noexcept
    {
        const auto r0 = this->r0(bs);
        const auto eps  = this->epsilon(bs);
        if(r <= r0)
        {
            return std::make_pair(-eps, 0.0);
        }
        const auto alp  = this->alpha(bs);
        const auto term = std::exp(-alp * (r - r0));

        return std::make_pair(eps * (term * term - 2 * term),
                              2 * alp * eps * term * (real_type(1) - term));
    }

  private:

    bool unit_converted_    = false;
    real_type K_BS_         = 6.0;
    real_type pi_over_K_BS_ = math::constants<real_type>::pi / 6.0;
    real_type alpha_BS_     = 3.0;

    std::array<real_type, 16> epsilon_BS_ = {{ // [kJ/mol]
        /* AA */ 14.39, /* AT */ 14.34, /* AG */ 13.25, /* AC */ 14.51,
        /* TA */ 10.37, /* TT */ 13.36, /* TG */ 10.34, /* TC */ 12.89,
        /* GA */ 14.81, /* GT */ 15.57, /* GG */ 14.93, /* GC */ 15.39,
        /* CA */ 11.42, /* CT */ 12.79, /* CG */ 10.52, /* CC */ 13.24
    }};
    std::array<real_type, 16> r0_BS_ = {{ // [angstrom]
        /* AA */ 3.716, /* AT */ 3.675, /* AG */ 3.827, /* AC */ 3.744,
        /* TA */ 4.238, /* TT */ 3.984, /* TG */ 4.416, /* TC */ 4.141,
        /* GA */ 3.576, /* GT */ 3.598, /* GG */ 3.664, /* GC */ 3.635,
        /* CA */ 4.038, /* CT */ 3.798, /* CG */ 4.208, /* CC */ 3.935
    }};
    std::array<real_type, 16> theta0_BS_ = {{ // [degree]
        /* AA */ 101.15, /* AT */ 85.94, /* AG */ 105.26, /* AC */ 89.00,
        /* TA */ 101.59, /* TT */ 89.50, /* TG */ 104.31, /* TC */ 91.28,
        /* GA */ 100.89, /* GT */ 84.83, /* GG */ 105.48, /* GC */ 88.28,
        /* CA */ 106.49, /* CT */ 93.31, /* CG */ 109.54, /* CC */ 95.46
    }};
};

} // mjolnir
#endif // MJOLNIR_POTENTIAL_LOCAL_3SPN2_BASE_STACKING_INTEARACTION_HPP
