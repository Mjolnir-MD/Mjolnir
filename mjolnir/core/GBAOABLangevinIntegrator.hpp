#ifndef MJOLNIR_CORE_G_BAOAB_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CORE_G_BAOAB_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>

namespace mjolnir
{

// g-BAOAB Langevin integrator developed by the following papers
// Leimkuhler B., Matthews C., Proc. R. Soc. A Math. Phys. Eng. Sci. (2016)
template<typename traitsT>
class GBAOABLangevinIntegrator
{
    public:
      using traits_type     = traitsT;
      using real_type       = typename traits_type::real_type;
      using coordinate_type = typename traits_type::coordinate_type;
      using system_type     = System<traitsT>;
      using forcefield_type = ForceField<traitsT>;
      using rng_type        = RandomNumberGenerator<traits_type>;

    public:

      GBAOABLangevinIntegrator(const real_type dt, std::vector<real_type>&& gamma)
          : dt_(dt), halfdt_(dt / 2), gammas_(std::move(gamma)),
            noise_coeff_ (gammas_.size())
      {}
      ~GBAOABLangevinIntegrator() = default;

      void initialize(system_type& sys, forcefield_type& ff, rng_type& rng);

      real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                     rng_type& rng);

    private:
      real_type dt_;
      real_type halfdt_;
      real_type temperature_;

      std::vector<real_type> gammas_;
      std::vector<real_type> noise_coeff_;
};

template<typename traitsT>
void GBAOABLangevinIntegrator<traitsT>::initialize(
        system_type& system, forcefield_type& ff, rng_type&)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    
    // calculate parameters for each particles
    this->update(system);

    // calculate force
    for(std::size_t i=0; i<system.size(); ++i)
    {
        system.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }
    ff.calc_force(system);
    return;
}

template<typename traitsT>
typename GBAOABLangevinIntegrator<traitsT>::real_type
GBAOABLangevinIntegrator<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
//    real_type largest_disp2(0.0);
//    for(std::size_t i=0; i<sys.size(); ++i)
//    {
//        auto& rm = sys.rmass(i);
//        auto& v  = sys.velocity(i);
//        auto& f  = sys.force(i);
//
//        v += this->halfdt_ * rm * f;
//    }
    return time + dt_;
}

} // mjolnir
#endif /* MJOLNIR_CORE_G_BAOAB_LANGEVIN_INTEGRATOR_HPP */
