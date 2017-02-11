#ifndef MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#define MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#include "Integrator.hpp"
#include "BoundaryCondition.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>

namespace mjolnir
{

template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class NVENewtonian : public Integrator<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef boundaryT boundary_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    NVENewtonian(const time_type dt, const std::size_t number_of_particles)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5),
          acceleration_(number_of_particles)
    {}
    ~NVENewtonian() override = default;

    void initialize(const ParticleContainer<traitsT>& pcon) override;
    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
                   ForceField<traitsT>& ff) override;

    time_type& delta_t()       override {return dt_;}
    time_type  delta_t() const override {return dt_;}

  private:
    time_type dt_;      //!< dt
    time_type halfdt_;  //!< dt/2
    time_type halfdt2_; //!< dt^2/2
    std::vector<coordinate_type> acceleration_;
};

template<typename traitsT, typename boundaryT>
void NVENewtonian<traitsT, boundaryT>::initialize(
        const ParticleContainer<traitsT>& pcon)
{
    for(auto iter = make_zip(pcon.cbegin(), acceleration_.begin());
            iter != make_zip(pcon.cend(), acceleration_.end()); ++iter)
    {
        *get<1>(iter) = get<0>(iter)->force / get<0>(iter)->mass;
    }
    return;
}

// at the initial step, acceleration_ must be initialized
template<typename traitsT, typename boundaryT>
typename NVENewtonian<traitsT, boundaryT>::time_type
NVENewtonian<traitsT, boundaryT>::step(const time_type time,
        ParticleContainer<traitsT>& pcon, ForceField<traitsT>& ff)
{
    // calc r(t+dt)
    for(auto iter = make_zip(pcon.begin(), acceleration_.cbegin());
            iter != make_zip(pcon.end(), acceleration_.cend()); ++iter)
    {
        get<0>(iter)->position = boundary_type::adjust_absolute(
            get<0>(iter)->position + dt_ * (get<0>(iter)->velocity) +
            halfdt2_ * (*get<1>(iter)));
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
        get<0>(iter)->force = coordinate_type(0., 0., 0.);
    }

    return time + dt_;
}

} // mjolnir
#endif /* MJOLNIR_VELOCITY_VERLET_INTEGRATOR */
