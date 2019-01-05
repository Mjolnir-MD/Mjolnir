#ifndef MJOLNIR_NVE_VELOCITY_VERLET_INTEGRATOR_HPP
#define MJOLNIR_NVE_VELOCITY_VERLET_INTEGRATOR_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>

namespace mjolnir
{

template<typename traitsT>
class VelocityVerletIntegrator
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef System<traitsT>     system_type;
    typedef ForceField<traitsT> forcefield_type;

  public:

    VelocityVerletIntegrator(const real_type dt) noexcept
        : dt_(dt), halfdt_(dt / 2)
    {}
    ~VelocityVerletIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff);

    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    real_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt / 2;
    }

    void update(const system_type& sys) const noexcept {/* do nothing */}

  private:
    real_type dt_;      //!< dt
    real_type halfdt_;  //!< dt/2
};

template<typename traitsT>
void VelocityVerletIntegrator<traitsT>::initialize(
        system_type& system, forcefield_type& ff)
{
    this->update(system);

    for(std::size_t i=0; i<system.size(); ++i)
    {
        system[i].force = coordinate_type(real_type(0.0), real_type(0.0), real_type(0.0));
    }
    system.largest_displacement() = 0;
    ff.calc_force(system);
    return;
}

template<typename traitsT>
typename VelocityVerletIntegrator<traitsT>::real_type
VelocityVerletIntegrator<traitsT>::step(
        const real_type time, system_type& system, forcefield_type& ff)
{
    real_type largest_disp2(0);
    for(std::size_t i=0; i<system.size(); ++i)
    {
        auto pv = system[i]; // particle_view that points i-th particle.

        pv.velocity += (halfdt_ * pv.rmass) * pv.force;

        const auto disp = dt_ * pv.velocity;

        pv.position = system.adjust_position(pv.position + disp);
        pv.force    = coordinate_type(real_type(0.0), real_type(0.0), real_type(0.0));

        largest_disp2 = std::max(largest_disp2, length_sq(disp));
    }
    system.largest_displacement() = std::sqrt(largest_disp2);

    // calc f(t+dt)
    ff.calc_force(system);

    // calc v(t+dt)
    for(std::size_t i=0; i<system.size(); ++i)
    {
        auto pv = system[i];
        pv.velocity += (halfdt_ * pv.rmass) * pv.force;
    }
    return time + dt_;
}

} // mjolnir
#endif // MJOLNIR_NVE_VELOCITY_VERLET_INTEGRATOR_HPP