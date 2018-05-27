#ifndef MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include "constants.hpp"
#include "RandomNumberGenerator.hpp"

namespace mjolnir
{

template<typename traitsT>
class UnderdampedLangevinStepper
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type>     system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef RandomNumberGenerator<traits_type>    rng_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::system_motion_remover system_motion_remover;

    struct parameter_set
    {
        real_type       gamma;
        real_type       r_mass;
        real_type       sqrt_gamma_over_mass;
        coordinate_type accel;
    };

  public:

    UnderdampedLangevinStepper(const real_type dt,
            std::vector<real_type>&& gamma, rng_type&& rng)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5),
          rng_(std::move(rng)), parameters_(gamma.size())
    {
        for(std::size_t i=0; i<gamma.size(); ++i)
        {
            parameters_[i].gamma = gamma[i];
        }
    }
    ~UnderdampedLangevinStepper() = default;

    void initialize(system_type& sys, forcefield_type& ff);
    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    real_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt * 0.5; halfdt2_ = dt * dt * 0.5;
    }

    void set_temperature(const real_type T)
    {
        this->noise_coef_ =
            std::sqrt(2 * physics<real_type>::kB * T / dt_);
        return;
    }

    rng_type&       random_number_generator()       noexcept {return rng_;}
    rng_type const& random_number_generator() const noexcept {return rng_;}

  private:

    coordinate_type gen_random_accel(const real_type sqrt_g_over_m)
    {
        const real_type coef = this->noise_coef_ * sqrt_g_over_m;
        return coordinate_type(this->rng_.gaussian(0., coef),
                               this->rng_.gaussian(0., coef),
                               this->rng_.gaussian(0., coef));
    }

  private:
    real_type dt_;
    real_type halfdt_;
    real_type halfdt2_;
    real_type noise_coef_;
    rng_type  rng_;

    std::vector<parameter_set> parameters_;
};

template<typename traitsT>
void UnderdampedLangevinStepper<traitsT>::initialize(
        system_type& system, forcefield_type& ff)
{
    system_motion_remover::invoke(system);

    this->noise_coef_ = std::sqrt(
        2 * physics<real_type>::kB * system.attribute("temperature") / dt_);

    for(std::size_t i=0; i<system.size(); ++i)
    {
        system[i].force = coordinate_type(0.0, 0.0, 0.0);
    }
    system.largest_displacement() = 0;
    ff.calc_force(system);

    for(std::size_t i=0; i<system.size(); ++i)
    {
        auto& p = parameters_[i];
        p.r_mass = 1.0 / system[i].mass;
        p.sqrt_gamma_over_mass = std::sqrt(p.gamma * p.r_mass);
        p.accel = system[i].force * p.r_mass +
                  this->gen_random_accel(p.sqrt_gamma_over_mass);
    }
    return;
}

template<typename traitsT>
typename UnderdampedLangevinStepper<traitsT>::real_type
UnderdampedLangevinStepper<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff)
{
    real_type largest_disp2(0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& param = this->parameters_[i];

        const real_type gamma_dt_over_2           = param.gamma * halfdt_;
        const real_type one_minus_gamma_dt_over_2 = 1. - gamma_dt_over_2;

        const auto disp = (dt_ * one_minus_gamma_dt_over_2) * (sys[i].velocity) +
                halfdt2_ * param.accel;

        sys[i].position = sys.adjust_position(sys[i].position + disp);

        sys[i].velocity *= one_minus_gamma_dt_over_2 *
            (one_minus_gamma_dt_over_2 * one_minus_gamma_dt_over_2 +
             gamma_dt_over_2);
        sys[i].velocity += (halfdt_ * one_minus_gamma_dt_over_2) * param.accel;

        sys[i].force = coordinate_type(0., 0., 0.);

        largest_disp2 = std::max(largest_disp2, length_sq(disp));
    }
    sys.largest_displacement() = std::sqrt(largest_disp2);

    // calc f(t+dt)
    ff.calc_force(sys);

    // calc a(t+dt) and v(t+dt), generate noise
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        auto& param = this->parameters_[i];
        param.accel = sys[i].force * param.r_mass +
                      this->gen_random_accel(param.sqrt_gamma_over_mass);

        sys[i].velocity += halfdt_ * (1. - param.gamma * halfdt_) * param.accel;
    }

    system_motion_remover::invoke(sys);

    return time + dt_;
}


} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
