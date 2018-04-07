#ifndef MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#define MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>

namespace mjolnir
{

template<typename traitsT>
class VelocityVerletStepper
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef System<traitsT>     system_type;
    typedef ForceField<traitsT> forcefield_type;
    typedef typename system_type::particle_type particle_type;

  public:

    VelocityVerletStepper(const real_type dt) noexcept
        : dt_(dt), halfdt_(dt * 0.5)
    {}
    ~VelocityVerletStepper() = default;

    void initialize(system_type& sys, forcefield_type& ff);

    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    real_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt * 0.5;
    }

  private:
    real_type dt_;      //!< dt
    real_type halfdt_;  //!< dt/2
};

template<typename traitsT>
void VelocityVerletStepper<traitsT>::initialize(
        system_type& system, forcefield_type& ff)
{
    for(std::size_t i=0; i<system.size(); ++i)
    {
        system[i].force = coordinate_type(0.0, 0.0, 0.0);
    }
    system.largest_displacement() = 0;
    ff.calc_force(system);
    return;
}

template<typename traitsT>
typename VelocityVerletStepper<traitsT>::real_type
VelocityVerletStepper<traitsT>::step(
        const real_type time, system_type& system, forcefield_type& ff)
{
    real_type largest_disp2(0);
    for(std::size_t i=0; i<system.size(); ++i)
    {
        auto& particle = system[i];

        particle.velocity += (halfdt_ / particle.mass) * particle.force;

        const auto disp = dt_ * particle.velocity;

        particle.position = system.adjust_position(particle.position + disp);
        particle.force    = coordinate_type(0.0, 0.0, 0.0);

        largest_disp2 = std::max(largest_disp2, length_sq(disp));
    }
    system.largest_displacement() = std::sqrt(largest_disp2);

    // calc f(t+dt)
    ff.calc_force(system);

    // calc v(t+dt)
    for(std::size_t i=0; i<system.size(); ++i)
    {
        auto& particle = system[i];
        particle.velocity += (halfdt_ / particle.mass) * particle.force;
    }

    remove_translation(system);
    remove_rotation(system);

    return time + dt_;
}

} // mjolnir
#endif /* MJOLNIR_VELOCITY_VERLET_INTEGRATOR */
