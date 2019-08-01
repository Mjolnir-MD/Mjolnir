#ifndef MJOLNIR_CORE_VELOCITY_VERLET_INTEGRATOR_HPP
#define MJOLNIR_CORE_VELOCITY_VERLET_INTEGRATOR_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>

namespace mjolnir
{

template<typename traitsT>
class VelocityVerletIntegrator
{
  public:
    using traits_type = traitsT;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traitsT>;
    using forcefield_type = ForceField<traitsT>;

  public:

    explicit VelocityVerletIntegrator(const real_type dt) noexcept
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

    void update(const system_type&) const noexcept {/* do nothing */}

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
        system.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }
    ff.calc_force(system);
    return;
}

template<typename traitsT>
typename VelocityVerletIntegrator<traitsT>::real_type
VelocityVerletIntegrator<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff)
{
    real_type largest_disp2(0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i) += (halfdt_ * sys.rmass(i)) * sys.force(i);

        const auto disp = dt_ * sys.velocity(i);

        sys.position(i) = sys.adjust_position(sys.position(i) + disp);
        sys.force(i)    = math::make_coordinate<coordinate_type>(0, 0, 0);

        largest_disp2 = std::max(largest_disp2, math::length_sq(disp));
    }

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff.update_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(t+dt)
    ff.calc_force(sys);

    // calc v(t+dt)
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i) += (halfdt_ * sys.rmass(i)) * sys.force(i);
    }
    return time + dt_;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif // MJOLNIR_NVE_VELOCITY_VERLET_INTEGRATOR_HPP
