#ifndef MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include "constants.hpp"
#include "RandomNumberGenerator.hpp"
#include <mjolnir/util/observer_ptr.hpp>
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
    typedef observer_ptr<RandomNumberGenerator<traits_type>> rng_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    UnderdampedLangevinStepper(const time_type dt,
        std::vector<real_type>&& gamma, RandomNumberGenerator<traitsT>& rng)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5), rng_(&rng),
          gamma_(std::move(gamma)), noise_(gamma_.size()), accel_(gamma_.size())
    {}
    ~UnderdampedLangevinStepper() = default;

    void initialize(const system_type& pcon);
    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    time_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt * 0.5; halfdt2_ = dt*dt*0.5;
    }

  private:

    coordinate_type gen_random_accel(const real_type m, const real_type g)
    {
        const real_type coef = this->noise_coef_ * std::sqrt(g / m);
        return coordinate_type(this->rng_->gaussian(0., coef),
                               this->rng_->gaussian(0., coef),
                               this->rng_->gaussian(0., coef));
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

template<typename traitsT, typename randomT>
void UnderdampedLangevinStepper<traitsT, randomT>::initialize(
        const system_type& system)
{
    this->noise_coef =
        std::sqrt(2 * physics<real_type>::kB * sys.temperature() / dt_);

    for(auto iter = make_zip(system.cbegin(), gamma_.cbegin(),
                             noise_.begin(),  accel_.begin());
            iter != make_zip(system.cend(),   gamma_.cend(),
                             noise_.end(),    accel_.end()); ++iter)
    {
        *get<3>(iter) = get<0>(iter)->force / get<0>(iter)->mass;
        *get<2>(iter) = gen_random_accel(get<0>(iter)->mass,
                            /*gamma = */*get<1>(iter));
    }

    return;

}

template<typename traitsT, typename randomT>
typename UnderdampedLangevinStepper<traitsT, randomT>::real_type
UnderdampedLangevinStepper<traitsT, randomT>::step(
        const real_type time, system_type& sys, forcefield_type& ff)
{
    this->noise_coef =
        std::sqrt(2 * physics<real_type>::kB * sys.temperature() / dt_);

    for(auto iter(make_zip(sys.begin(),     gamma_.cbegin(),
                           noise_.cbegin(), accel_.cbegin())),
             end_(make_zip(sys.end(),       gamma_.cend(),
                           noise_.cend(),   accel_.cend())),
            iter != end_; ++iter)
    {
        const real_type hgdt   = (*get<1>(iter) * halfdt_);
        const real_type o_hgdt = 1. - hgdt;

        const coordinate_type noisy_force = ((*get<2>(iter)) + (*get<3>(iter)) / m);

        get<0>(iter)->position = boundary_type::adjust_absolute(
                get<0>(iter)->position +
                (dt_ * o_hgdt) * (get<0>(iter)->velocity) +
                halfdt2_ * noisy_force);

        get<0>(iter)->velocity *= o_hgdt * (o_hgdt * o_hgdt + hgdt);
        get<0>(iter)->velocity += (halfdt_ * o_hgdt) * noisy_force;
    }

    // calc f(t+dt)
    ff.calc_force(pcon);

    // calc a(t+dt) and v(t+dt), generate noise
    for(auto iter(make_zip(sys.begin(),    gamma_.cbegin(),
                           noise_.begin(), accel_.begin()));
             end_(make_zip(sys.end(),      gamma_.cend(),
                           noise_.end(),   accel_.end()));
            iter != end_; ++iter)
    {
        // consider cash 1/m
        const coordinate_type acc = get<0>(iter)->force / (get<0>(iter)->mass);
        *get<3>(iter) = acc;

        const coordinate_type noise =
            this->gen_random_accel(get<0>(iter)->mass, *get<1>(iter));
        *get<2>(iter) = noise;

        const real_type gm = (*get<1>(iter));
        get<0>(iter)->velocity +=
            halfdt_ * (1. - gm * halfdt_) * (acc + noise);
        get<0>(iter)->force = coordinate_type(0., 0., 0.);
    }

    return time + dt_;
}


} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
