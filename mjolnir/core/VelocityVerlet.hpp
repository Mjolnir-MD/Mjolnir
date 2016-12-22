#ifndef MJOLNIR_VELOCITY_VERLET_INTEGRATOR
#define MJOLNIR_VELOCITY_VERLET_INTEGRATOR
#include "Integrator.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>

namespace mjolnir
{

template<typename traitsT>
class VelocityVerlet : public Integrator<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    VelocityVerlet(const time_type dt, const std::size_t number_of_particles)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5),
          acceleration_(number_of_particles)
    {}
    ~VelocityVerlet() override = default;

    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
                   ForceField<traitsT>& ff) override;

  private:
    time_type dt_;      //!< dt
    time_type halfdt_;  //!< dt/2
    time_type halfdt2_; //!< dt^2/2
    std::vector<coordinate_type> acceleration_;
};

// at the initial step, acceleration_ must be initialized
template<typename traitsT>
typename VelocityVerlet<traitsT>::time_type
VelocityVerlet<traitsT>::step(const time_type time,
        ParticleContainer<traitsT>& pcon, ForceField<traitsT>& ff)
{
    // calc r(t+dt)
    for(auto iter = make_zip(pcon.begin(), acceleration_.begin());
            iter != make_zip(pcon.end(), acceleration_.end()); ++iter)
    {
        get<0>(iter)->position += dt_ * (get<0>(iter)->velocity) +
                                  halfdt2_ * (*get<1>(iter));
        get<0>(iter)->velocity += halfdt_ * (*get<1>(iter));
    }

    // calc f(t+dt)
    ff.calc_force(pcon);

    // calc a(t+dt) and v(t+dt)
    for(auto iter = make_zip(pcon.begin(), acceleration_.begin());
            iter != make_zip(pcon.end(), acceleration_.end()); ++iter)
    {
        // consider cash 1/m
        const coordinate_type acc = get<0>(iter)->force / (get<0>(iter)->mass);
        *get<1>(iter) = acc;
        get<0>(iter)->velocity += halfdt_ * acc;
    }

    return time + dt_;
}

} // mjolnir
#endif /* MJOLNIR_VELOCITY_VERLET_INTEGRATOR */
