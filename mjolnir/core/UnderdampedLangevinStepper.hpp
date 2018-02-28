#ifndef MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include "constants.hpp"
#include "RandomNumberGenerator.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>

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

  public:

    UnderdampedLangevinStepper(const real_type dt,
            std::vector<real_type>&& gamma, rng_type&& rng)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5), rng_(std::move(rng)),
          gamma_(std::move(gamma)), noise_(gamma_.size()), accel_(gamma_.size())
    {}
    ~UnderdampedLangevinStepper() = default;

    void initialize(const system_type& sys);
    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    real_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt * 0.5; halfdt2_ = dt*dt*0.5;
    }

    rng_type&       random_number_generator()       noexcept {return rng_;}
    rng_type const& random_number_generator() const noexcept {return rng_;}

  private:

    coordinate_type gen_random_accel(const real_type m, const real_type g)
    {
        const real_type coef = this->noise_coef_ * std::sqrt(g / m);
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
    std::vector<real_type>       gamma_;
    std::vector<coordinate_type> noise_;
    std::vector<coordinate_type> accel_;
};

template<typename traitsT>
void UnderdampedLangevinStepper<traitsT>::initialize(
        const system_type& system)
{
    this->noise_coef_ =
        std::sqrt(2 * physics<real_type>::kB * system.temperature() / dt_);

    for(std::size_t i=0; i<system.size(); ++i)
    {
        accel_[i] = system[i].force / system[i].mass;
        noise_[i] = this->gen_random_accel(system[i].mass, gamma_[i]);
    }
    return;
}

template<typename traitsT>
typename UnderdampedLangevinStepper<traitsT>::real_type
UnderdampedLangevinStepper<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff)
{
    this->noise_coef_ =
        std::sqrt(2 * physics<real_type>::kB * sys.temperature() / dt_);

    real_type max_speed2(0.);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        // max of v(t)
        max_speed2 = std::max(max_speed2, length_sq(sys[i].velocity));

        const real_type hgdt   = gamma_[i] * halfdt_;
        const real_type o_hgdt = 1. - hgdt;

        const coordinate_type noisy_force = noise_[i] + accel_[i];

        sys[i].position = sys.adjust_position(sys[i].position +
                (dt_ * o_hgdt) * (sys[i].velocity) +
                halfdt2_ * noisy_force);

        sys[i].velocity *= o_hgdt * (o_hgdt * o_hgdt + hgdt);
        sys[i].velocity += (halfdt_ * o_hgdt) * noisy_force;
        // here, v comes v(t+h/2)
    }

    sys.max_speed() = std::sqrt(max_speed2);

    // calc f(t+dt)
    ff.calc_force(sys);

    // calc a(t+dt) and v(t+dt), generate noise
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        // consider cash 1/m
        const coordinate_type acc = sys[i].force / sys[i].mass;
        accel_[i] = acc;

        const coordinate_type noise =
            this->gen_random_accel(sys[i].mass, gamma_[i]);
        noise_[i] = noise;

        sys[i].velocity += halfdt_ * (1. - gamma_[i] * halfdt_) * (acc + noise);
        sys[i].force = coordinate_type(0., 0., 0.);
    }

    return time + dt_;
}


} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
