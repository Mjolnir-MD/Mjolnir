#ifndef MJOLNIR_CORE_VELOCITY_VERLET_INTEGRATOR_HPP
#define MJOLNIR_CORE_VELOCITY_VERLET_INTEGRATOR_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/util/logger.hpp>

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
    using matrix33_type   = typename traits_type::matrix33_type;
    using system_type     = System<traitsT>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traitsT>>;
    using rng_type        = RandomNumberGenerator<traitsT>;
    using remover_type    = SystemMotionRemover<traits_type>;

  public:

    explicit VelocityVerletIntegrator(
            const real_type dt, remover_type&& remover) noexcept
        : dt_(dt), halfdt_(dt / 2), remover_(std::move(remover))
    {}
    ~VelocityVerletIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type& rng);

    real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                   rng_type& rng);

    real_type delta_t() const noexcept {return dt_;}

    void update(const system_type&) noexcept {/* do nothing */}
    void update(const system_type&, const real_type newdt) noexcept
    {
        this->dt_     = newdt;
        this->halfdt_ = newdt * 0.5;
        return ;
    }

  private:
    real_type dt_;      //!< dt
    real_type halfdt_;  //!< dt/2

    remover_type remover_;
};

template<typename traitsT>
void VelocityVerletIntegrator<traitsT>::initialize(
        system_type& system, forcefield_type& ff, rng_type&)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if(!ff->constraint().empty())
    {
        MJOLNIR_LOG_WARN(
            "Velocity verlet integrator does not support constraint forcefield."
            " [[forcefields.constraint]] will be ignored.");
    }

    this->update(system);

    // if loaded from MsgPack, we can skip it.
    if( ! system.force_initialized())
    {
        for(std::size_t i=0; i<system.size(); ++i)
        {
            system.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        for(auto& kv : system.variables())
        {
            auto& var = kv.second;
            var.update(var.x(), var.v(), real_type(0));
        }
        system.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);
        ff->calc_force(system);
    }
    return;
}

template<typename traitsT>
typename VelocityVerletIntegrator<traitsT>::real_type
VelocityVerletIntegrator<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff, rng_type&)
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
    for(auto& kv : sys.variables())
    {
        auto& var = kv.second;
        auto next_x = var.x();
        auto next_v = var.v();

        next_v += halfdt_ * var.f() / var.m();
        next_x += dt_ * next_v;
        var.update(next_x, next_v, real_type(0));
    }
    sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff->reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(t+dt)
    ff->calc_force(sys);

    // calc v(t+dt)
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i) += (halfdt_ * sys.rmass(i)) * sys.force(i);
    }
    for(auto& kv : sys.variables())
    {
        auto& var = kv.second;
        var.update(var.x(), var.v() + halfdt_ * var.f() / var.m(), var.f());
    }

    // remove net rotation/translation
    remover_.remove(sys);

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
