#ifndef MJOLNIR_CORE_G_BAOAB_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CORE_G_BAOAB_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// g-BAOAB Langevin integrator developed by the following papers
// Leimkuhler B., Matthews C., Proc. R. Soc. A Math. Phys. Eng. Sci. (2016)
template<typename traitsT>
class gBAOABLangevinIntegrator
{
  public:
    using traits_type                = traitsT;
    using real_type                  = typename traits_type::real_type;
    using coordinate_type            = typename traits_type::coordinate_type;
    using indices_type               = std::array<std::size_t, 2>;
    using system_type                = System<traitsT>;
    using forcefield_type            = ForceField<traitsT>;
    using constraint_forcefield_type = typename forcefield_type::constraint_forcefield_type;
    using constraint_type            = std::pair<indices_type, real_type>;
    using constraints_type           = std::vector<constraint_type>;
    using rng_type                   = RandomNumberGenerator<traits_type>;
    using remover_type               = SystemMotionRemover<traits_type>;

  public:

    gBAOABLangevinIntegrator(const real_type dt, std::vector<real_type>&& gamma,
                             remover_type&& remover)
        : dt_(dt), halfdt_(dt / 2), gammas_(std::move(gamma)),
          exp_gamma_dt_(gammas_.size()), noise_coeff_ (gammas_.size()),
          remover_(std::move(remover)),
          old_position_(gammas_.size()),
          old_pos_rattle_(gammas_.size())
    {}
    ~gBAOABLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type& rng);

    real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                   rng_type& rng);

    void update(const system_type& sys)
    {
        if(!sys.has_attribute("temperature"))
        {
            throw std::out_of_range("mjolnir::g-BAOABLangevinIntegrator: "
                "Langevin Integrator requires reference temperature, but "
                "`temperature` is not found in `system.attribute`.");
        }
        this->temperature_ = sys.attribute("temperature");
        this->reset_parameter(sys);
        return;
    }

    real_type delta_t() const noexcept {return dt_;}
    std::vector<real_type> const& parameters() const noexcept {return gammas_;}

  private:

    void reset_parameter(const system_type& sys) noexcept
    {
        const auto kBT = physics::constants<real_type>::kB() * this->temperature_;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto gamma     = this->gammas_.at(i);
            const auto gamma_dt  = -1 * gamma * this->dt_;
            this->exp_gamma_dt_.at(i) = std::exp(gamma_dt);
            this->noise_coeff_ .at(i) = std::sqrt(
                    kBT * (1 - std::exp(2 * gamma_dt)) * sys.rmass(i));
        }
        return;
    };

    coordinate_type gen_R(rng_type& rng) noexcept
    {
        const auto x = rng.gaussian();
        const auto y = rng.gaussian();
        const auto z = rng.gaussian();
        return math::make_coordinate<coordinate_type>(x, y, z);
    }

    void correct_coordinate(system_type& sys, const forcefield_type& ff)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        std::size_t rattle_step = 0;
        const constraint_forcefield_type& constraint_ff = ff.constraint();
        const constraints_type&           constraints = constraint_ff.constraints();

        while(rattle_step<constraint_ff.max_iteration())
        {
            bool corrected = false;
            for(std::size_t i=0; i<constraints.size(); ++i)
            {
                const auto& indices_v0 = constraints[i];
                const auto& indices    = indices_v0.first;
                auto& p1 = sys.position(indices[0]);
                auto& p2 = sys.position(indices[1]);
                auto& v1 = sys.velocity(indices[0]);
                auto& v2 = sys.velocity(indices[1]);

                const auto      dp = sys.adjust_direction(p2 - p1);
                const real_type dp2 = math::length_sq(dp);
                const real_type missmatch2 = square_v0s_[i] - dp2;

                if(correction_tolerance_ < std::abs(missmatch2))
                {
                    auto& op1 = this->old_pos_rattle_[indices[0]];
                    auto& op2 = this->old_pos_rattle_[indices[1]];

                    const auto old_dp = sys.adjust_direction(op2 - op1);
                    const real_type dot_old_new_dp = math::dot_product(old_dp, dp);
                    const real_type lambda =
                        0.5 * missmatch2 * reduced_mass_[i] / dot_old_new_dp;
                    const coordinate_type correction_force = lambda * old_dp;
                    const auto& rm1 = sys.rmass(indices[0]);
                    const auto& rm2 = sys.rmass(indices[1]);
                    const coordinate_type correction_vec1  = correction_force * rm1;
                    const coordinate_type correction_vec2  = correction_force * rm2;
                    p1 -= correction_vec1;
                    p2 += correction_vec2;
                    v1 -= correction_vec1 * r_dt_in_correction_;
                    v2 += correction_vec2 * r_dt_in_correction_;

                    corrected = true;
                }
            }

            if(!corrected) {break;}

            ++rattle_step;
        }

        if(rattle_step >= constraint_ff.max_iteration())
        {
            MJOLNIR_LOG_WARN("coordinate rattle iteration number exceeds rattle max iteration");
        }

        return;
    }

    void correct_velocity(system_type& sys, const forcefield_type& ff)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        const constraint_forcefield_type& constraint_ff = ff.constraint();
        const constraints_type&           constraints   = constraint_ff.constraints();

        std::size_t rattle_step = 0;
        while(rattle_step<constraint_ff.max_iteration())
        {
            bool corrected = false;
            for(std::size_t i=0; i<constraints.size(); ++i)
            {
                const auto& constraint       = constraints[i];
                const auto& indices          = constraint.first;
                auto& p1  = sys.position(indices[0]);
                auto& p2  = sys.position(indices[1]);
                auto& v1  = sys.velocity(indices[0]);
                auto& v2  = sys.velocity(indices[1]);

                const auto pos_diff      = sys.adjust_direction(p2 - p1);
                const auto vel_diff      = v2 - v1;
                const real_type dot_pdvd = math::dot_product(pos_diff, vel_diff);
                const real_type lambda   = dot_pdvd * reduced_mass_[i] * r_square_v0s_[i];

                if(correction_tolerance_ < std::abs(lambda))
                {
                    const auto& rm1 = sys.rmass(indices[0]);
                    const auto& rm2 = sys.rmass(indices[1]);

                    const auto correction_vec = lambda * pos_diff;
                    v1 += correction_vec * rm1;
                    v2 -= correction_vec * rm2;
                    corrected = true;
                }
            }

            if(!corrected) {break;}

            ++rattle_step;
        }

        if(constraint_ff.max_iteration() <= rattle_step)
        {
            MJOLNIR_LOG_WARN("velocity rattle iteration number exceeds rattle max iteration.");
        }
        return;
   }

  private:
    real_type   dt_;
    real_type   halfdt_;
    std::vector<real_type>       gammas_;
    std::vector<real_type>       exp_gamma_dt_;
    std::vector<real_type>       noise_coeff_;
    std::vector<real_type>       square_v0s_;
    std::vector<real_type>       r_square_v0s_;
    std::vector<real_type>       reduced_mass_;
    remover_type remover_;

    real_type   temperature_;
    std::size_t correction_max_iter_;
    real_type   correction_tolerance_;
    real_type   correction_tolerance_dt_;
    real_type   correction_tolerance_dt_itr_;
    real_type   dt_in_correction_;
    real_type   r_dt_in_correction_;

    static constexpr std::size_t correction_iter_num_ = 1;

    std::vector<coordinate_type> old_position_;
    std::vector<coordinate_type> old_pos_rattle_;

};

template<typename traitsT>
void gBAOABLangevinIntegrator<traitsT>::initialize(
        system_type& system, forcefield_type& ff, rng_type&)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
  
    const constraint_forcefield_type& constraint_ff = ff.constraint();
    const constraints_type&           constraints   = constraint_ff.constraints();
    const real_type                   tolerance     = constraint_ff.tolerance();
    // calculate parameters for constraint
    this->correction_tolerance_        = tolerance;
    this->correction_tolerance_dt_     = tolerance / dt_;
    this->correction_tolerance_dt_itr_ = tolerance * 2. * correction_iter_num_ / dt_;
    this->dt_in_correction_            = dt_ * 0.5 / correction_iter_num_;
    this->r_dt_in_correction_          = 1. / dt_in_correction_;


    // calculate parameters for each particles
    this->update(system);

    // calculate force
    for(std::size_t i=0; i<system.size(); ++i)
    {
        system.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }
    ff.calc_force(system);

    // buffering old position
    for(std::size_t i=0; i<system.size(); ++i)
    {
        old_position_[i] = system.position(i);
    }

    MJOLNIR_LOG_INFO("hoge");
    square_v0s_.resize(constraints.size());
    r_square_v0s_.resize(constraints.size());
    reduced_mass_.resize(constraints.size());
    for(std::size_t i=0; i<constraints.size(); ++i)
    {
        // calculate square v0
        square_v0s_[i] = std::pow(constraints[i].second, 2);
        r_square_v0s_[i] = 1.0 / square_v0s_[i];

        // calculate inverse reduced mass
        const auto& constraint = constraints[i];
        std::size_t first_idx  = constraint.first[0];
        std::size_t second_idx = constraint.first[1];
        reduced_mass_[i] = 1.0 / (system.rmass(first_idx) + system.rmass(second_idx));
    }
    MJOLNIR_LOG_INFO("huga");

    return;
}

template<typename traitsT>
typename gBAOABLangevinIntegrator<traitsT>::real_type
gBAOABLangevinIntegrator<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    // B step
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i)  += this->halfdt_ * sys.rmass(i) * sys.force(i);
    }
    correct_velocity(sys, ff);

    // A step
    for(std::size_t correction_step=0; correction_step<correction_iter_num_; ++correction_step)
    {
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->old_pos_rattle_[i] = sys.position(i);
            sys.position(i) += this->dt_in_correction_ * sys.velocity(i);
        }
        correct_coordinate(sys, ff);
        correct_velocity(sys, ff);
    }

    // O step
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i) *= this->exp_gamma_dt_[i]; // *= exp(- gamma dt)
        sys.velocity(i) += this->noise_coeff_[i] * this->gen_R(rng);
    }
    correct_velocity(sys, ff);

    // A step
    for(std::size_t correction_step=0; correction_step<correction_iter_num_; ++correction_step)
    {
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            this->old_pos_rattle_[i] = sys.position(i);
            sys.position(i) += this->dt_in_correction_ * sys.velocity(i);
        }
        correct_coordinate(sys, ff);
        correct_velocity(sys, ff);
    }
    // update neighbor list; reduce margin, reconstruct the list if needed;
    real_type largest_disp2(0.0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        coordinate_type displacement = sys.position(i) - this->old_position_[i];
        largest_disp2 = std::max(largest_disp2, math::length_sq(displacement));
        sys.position(i) = sys.adjust_position(sys.position(i));

        // reset force
        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }

    ff.reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // B step
    ff.calc_force(sys);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i) += this->halfdt_ * sys.rmass(i) * sys.force(i);
    }
    correct_velocity(sys, ff);

    remover_.remove(sys);

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        this->old_position_[i] = sys.position(i);
    }

    return time + dt_;
}

} // mjolnir
#endif /* MJOLNIR_CORE_G_BAOAB_LANGEVIN_INTEGRATOR_HPP */
