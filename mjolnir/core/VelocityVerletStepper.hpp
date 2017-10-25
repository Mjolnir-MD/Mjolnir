#ifndef MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#define MJOLNIR_NVE_NEWTONIAN_INTEGRATOR
#include <mjolnir/core/System.hpp>
// #include <mjolnir/core/Rescaler.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/util/make_zip.hpp>

namespace mjolnir
{

template<typename traitsT>
struct NVE
{
    // currently, empty
};

template<typename traitsT, typename rescalingT = NVE<traitsT>>
class VelocityVerletStepper
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef rescalingT          rescaler_type;
    typedef System<traitsT>     system_type;
    typedef ForceField<traitsT> forcefield_type;
    typedef typename system_type::particle_type particle_type;

  public:

    VelocityVerletStepper(const real_type dt) noexcept
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5)
    {}
    ~VelocityVerletStepper() = default;

    void initialize(const system_type& sys);

    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    real_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt * 0.5; halfdt2_ = halfdt_ * dt_;
    }

  private:
    real_type dt_;      //!< dt
    real_type halfdt_;  //!< dt/2
    real_type halfdt2_; //!< dt^2/2
    std::vector<coordinate_type> acceleration_;
};

template<typename traitsT, typename rescalingT>
void VelocityVerletStepper<traitsT, rescalingT>::initialize(const system_type& sys)
{
    acceleration_.resize(sys.size());
    std::transform(sys.cbegin(), sys.cend(), acceleration_.begin(),
            [](const particle_type& p){return p.force / p.mass;});
    return;
}

template<typename traitsT, typename rescalingT>
typename VelocityVerletStepper<traitsT, rescalingT>::real_type
VelocityVerletStepper<traitsT, rescalingT>::step(
        const real_type time, system_type& system, forcefield_type& ff)
{
    real_type max_speed2(0.);
    // calc r(t+dt)
    for(auto iter = make_zip(system.begin(), acceleration_.cbegin());
            iter != make_zip(system.end(),   acceleration_.cend()); ++iter)
    {
        const auto& particle = *(get<0>(iter));
        const auto& acc      = *(get<1>(iter));

        max_speed2 = std::max(max_speed2, length_sq(particle.velocity));

        get<0>(iter)->position = system.adjust_position(
            particle.position + dt_ * (particle.velocity) + halfdt2_ * acc);

        get<0>(iter)->velocity += halfdt_ * acc;
    }
    system.max_speed() = std::sqrt(max_speed2);

    // calc f(t+dt)
    ff.calc_force(system);

    // calc a(t+dt) and v(t+dt)
    for(auto iter = make_zip(system.begin(), acceleration_.begin());
            iter != make_zip(system.end(),   acceleration_.end()); ++iter)
    {
        const auto& particle = *(get<0>(iter));

        const coordinate_type acc = particle.force / (particle.mass);
       *get<1>(iter) = acc;
        get<0>(iter)->velocity += halfdt_ * acc;
        get<0>(iter)->force     = coordinate_type(0., 0., 0.);
    }

    remove_translation(system);
    remove_rotation(system);

    return time + dt_;
}



} // mjolnir
#endif /* MJOLNIR_VELOCITY_VERLET_INTEGRATOR */
