#ifndef MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include "Integrator.hpp"
#include "BoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <memory>

namespace mjolnir
{

template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class UnderdampedLangevin : public Integrator<traitsT>
{
  public:
    typedef traitsT   traits_type;
    typedef boundaryT boundary_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    UnderdampedLangevin(const time_type dt, const std::size_t number_of_particles,
            const real_type temperature, const real_type kB,
            std::vector<real_type>&& friction_constant,
            const std::shared_ptr<RandomNumberGenerator<traits_type>>& rng)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5), rng_(rng),
          temperature_(temperature), kB_(kB),
          gamma_(std::forward<std::vector<real_type>>(friction_constant)),
          noise_(number_of_particles), acceleration_(number_of_particles)
    {}
    ~UnderdampedLangevin() override = default;

    void initialize(const ParticleContainer<traitsT>& pcon) override;

    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
                   ForceField<traitsT>& ff) override;

    real_type  T() const {return temperature_;}
    real_type& T()       {return temperature_;}

    time_type  delta_t() const override {return dt_;}
    time_type& delta_t()       override {return dt_;}

  private:
    time_type dt_;      //!< dt
    time_type halfdt_;  //!< dt/2
    time_type halfdt2_; //!< dt^2/2
    real_type kB_;
    real_type temperature_;
    std::shared_ptr<RandomNumberGenerator<traits_type>> rng_;
    std::vector<real_type> gamma_;
    std::vector<coordinate_type> noise_;
    std::vector<coordinate_type> acceleration_;
};

template<typename traitsT, typename boundaryT>
void UnderdampedLangevin<traitsT, boundaryT>::initialize(
        const ParticleContainer<traitsT>& pcon)
{
    for(auto iter = make_zip(pcon.cbegin(),  gamma_.cbegin(),
                             noise_.begin(), acceleration_.begin());
            iter != make_zip(pcon.cend(),    gamma_.cend(),
                             noise_.end(),   acceleration_.end()); ++iter)
    {
        // set a = f/m
        *get<3>(iter) = get<0>(iter)->force / get<0>(iter)->mass;

        // set random force
        *get<2>(iter) = rng_->underdamped_langevin(
                get<0>(iter)->mass, *get<1>(iter), dt_, temperature_, kB_);
    }

    return;
}

// at the initial step, acceleration_ and noise_ must be initialized
template<typename traitsT, typename boundaryT>
typename UnderdampedLangevin<traitsT, boundaryT>::time_type
UnderdampedLangevin<traitsT, boundaryT>::step(const time_type time,
        ParticleContainer<traitsT>& pcon, ForceField<traitsT>& ff)
{
    // calc r(t+dt)
    for(auto iter = make_zip(pcon.begin(),    gamma_.cbegin(),
                             noise_.cbegin(), acceleration_.cbegin());
            iter != make_zip(pcon.end(),      gamma_.cend(),
                             noise_.cend(),   acceleration_.cend()); ++iter)
    {
        const real_type hgdt   = (*get<1>(iter) * halfdt_);
        const real_type o_hgdt = 1. - hgdt;

        const coordinate_type noisy_force = ((*get<2>(iter)) + (*get<3>(iter)));

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
    for(auto iter = make_zip(pcon.begin(),   gamma_.cbegin(),
                             noise_.begin(), acceleration_.begin());
            iter != make_zip(pcon.end(),     gamma_.cend(),
                             noise_.end(),   acceleration_.end()); ++iter)
    {
        // consider cash 1/m
        const coordinate_type acc = get<0>(iter)->force / (get<0>(iter)->mass);
        *get<3>(iter) = acc;

        const coordinate_type noise = rng_->underdamped_langevin(
                get<0>(iter)->mass, *get<1>(iter), dt_, temperature_, kB_);
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
