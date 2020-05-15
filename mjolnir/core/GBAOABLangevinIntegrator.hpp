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
      using traits_type      = traitsT;
      using real_type        = typename traits_type::real_type;
      using coordinate_type  = typename traits_type::coordinate_type;
      using indices_type     = std::array<std::size_t, 2>;
      using constraints_type = std::vector<std::pair<indices_type, real_type>>;
      using system_type      = System<traitsT>;
      using forcefield_type  = ForceField<traitsT>;
      using rng_type         = RandomNumberGenerator<traits_type>;
      using remover_type     = SystemMotionRemover<traits_type>;

    public:

      GBAOABLangevinIntegrator(const real_type dt, std::vector<real_type>&& gamma,
                               constraints_type&& constraint, remover_type&& remover)
          : dt_(dt), halfdt_(dt / 2), gammas_(std::move(gamma)),
            exp_gamma_dt_(gammas_.size()), noise_coeff_ (gammas_.size()), 
            remover_(std::move(remover))
      {}
      ~GBAOABLangevinIntegrator() = default;

      void initialize(system_type& sys, forcefield_type& ff, rng_type& rng);

      real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                     rng_type& rng);

      void update(const system_type& sys)
      {
          if(!sys.has_attribute("temperature"))
          {
              throw std::out_of_range("mjolnir::g-BAOABLangevinIntegrator: "
                  "Langevin Integrator requires reference temperature, but "
                  "`temperature` is not found in `system.attribute`.");
          }
          this->temperature_ = sys.attribute("temperature");
          this->reset_parameter(sys);
          return;
      }

      real_type delta_t() const noexcept {return dt_;}
      std::vector<real_type> const& parameters() const noexcept {return gammas_;}

    private:

      void reset_parameter(const system_type& sys) noexcept
      {
          const auto kBT = physics::constants<real_type>::kB() * this->temperature_;
          for(std::size_t i=0; i<sys.size(); ++i)
          {
              const auto gamma     = this->gammas_.at(i);
              const auto gamma_dt  = -1 * gamma * this->dt_;
              this->exp_gamma_dt_.at(i) = std::exp(gamma_dt);
              this->noise_coeff_ .at(i) = std::sqrt(
                      kBT * (1 - std::exp(2 * gamma_dt)) * sys.rmass(i));
          }
          return;
      };

      coordinate_type gen_R(rng_type& rng) noexcept
      {
          const auto x = rng.gaussian();
          const auto y = rng.gaussian();
          const auto z = rng.gaussian();
          return math::make_coordinate<coordinate_type>(x, y, z);
      }

    private:
      real_type dt_;
      real_type halfdt_;
      real_type temperature_;

      std::vector<real_type> gammas_;
      std::vector<real_type> exp_gamma_dt_;
      std::vector<real_type> noise_coeff_;
      std::vector<coordinate_type> old_position_;

      remover_type remover_;
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

    // buffering old position
    old_position_.reserve(system.size());
    for(std::size_t i=0; i<system.size(); ++i)
    {
        old_position_.push_back(system.position(i));
    }

    return;
}

template<typename traitsT>
typename GBAOABLangevinIntegrator<traitsT>::real_type
GBAOABLangevinIntegrator<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    // B step
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i)  += this->halfdt_ * sys.rmass(i) * sys.force(i);
    }
    // velocity correction step
    // TODO
    // A step
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.position(i) += this->halfdt_ * sys.velocity(i);
    }
    // coordinate correction step
    // TODO
    // velocity correction step
    // TODO
    // O step
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i) *= this->exp_gamma_dt_[i]; // *= exp(- gamma dt)
        sys.velocity(i) += this->noise_coeff_[i] * this->gen_R(rng);
    }
    // velocity correction step
    // TODO
    // A step
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.position(i) += this->halfdt_ * sys.velocity(i);
    }
    // coordinate correction step
    // TODO
    // velocity correction step
    // TODO
    // update neighbor list; reduce margin, reconstruct the list if needed;
    real_type largest_disp2(0.0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        coordinate_type displacement = sys.position(i) - this->old_position_[i];
        largest_disp2 = std::max(largest_disp2, math::length_sq(displacement));
        sys.position(i) = sys.adjust_position(sys.position(i));

        // reset force
        sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }

    ff.reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // B step
    ff.calc_force(sys);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        sys.velocity(i) += this->halfdt_ * sys.rmass(i) * sys.force(i);
    }
    // velocity correction step
    // TODO

    remover_.remove(sys);

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        this->old_position_[i] = sys.position(i);
    }

    return time + dt_;
}

} // mjolnir
#endif /* MJOLNIR_CORE_G_BAOAB_LANGEVIN_INTEGRATOR_HPP */
