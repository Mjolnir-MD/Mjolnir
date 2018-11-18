#ifndef MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/constants.hpp>

namespace mjolnir
{

// Simple underdamped Langevin integrator.
// The implementation is based on the article written by J. D. Honeycutt and
// D. Thirumalai (1992) Biopolymers and also the article written by Z. Guo and
// D. Thirumalai (1995) Biopolymers. Same algorithm used as default Langevin
// integrator in CafeMol that is developed by H. Kenzaki et al., (2011) JCTC.
template<typename traitsT>
class UnderdampedLangevinStepper
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type>     system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef RandomNumberGenerator<traits_type>    rng_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    struct parameter_set
    {
        real_type       gamma;
        real_type       sqrt_gamma_over_mass;
        coordinate_type accel;
    };

  public:

    UnderdampedLangevinStepper(const real_type dt,
            std::vector<real_type>&& gamma, rng_type&& rng)
        : dt_(dt), halfdt_(dt / 2), halfdt2_(dt * dt / 2),
          rng_(std::move(rng)), parameters_(gamma.size())
    {
        for(std::size_t i=0; i<gamma.size(); ++i)
        {
            parameters_[i].gamma = gamma[i];
        }
    }
    ~UnderdampedLangevinStepper() = default;

    void initialize(system_type& sys, forcefield_type& ff);
    real_type step(const real_type time, system_type& sys, forcefield_type& ff);

    real_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt / 2; halfdt2_ = dt * dt / 2;
    }

    void update(const system_type& sys)
    {
        if(!sys.has_attribute("temperature"))
        {
            throw std::out_of_range("mjolnir::UnderdampedLangevinStepper: "
                "Langevin Integrator requires reference temperature, but "
                "`temperature` is not found in `system.attribute`.");
        }

        this->temperature_ = sys.attribute("temperature");
        this->noise_coef_  = std::sqrt(2 * physics::constants<real_type>::kB() *
                                       this->temperature_ / dt_);
    }

  private:

    coordinate_type gen_gaussian_vec(const real_type coef)
    {
        return coordinate_type(this->rng_.gaussian(0.0, coef),
                               this->rng_.gaussian(0.0, coef),
                               this->rng_.gaussian(0.0, coef));
    }

  private:
    real_type dt_;
    real_type halfdt_;
    real_type halfdt2_;
    real_type temperature_;
    real_type noise_coef_;
    rng_type  rng_;

    std::vector<parameter_set> parameters_;
};

template<typename traitsT>
void UnderdampedLangevinStepper<traitsT>::initialize(
        system_type& system, forcefield_type& ff)
{
    this->update(system);

    for(std::size_t i=0; i<system.size(); ++i)
    {
        system[i].force = coordinate_type(real_type(0.0), real_type(0.0), real_type(0.0));
    }
    system.largest_displacement() = real_type(0.0);
    ff.calc_force(system);

    for(std::size_t i=0; i<system.size(); ++i)
    {
        auto pv = system[i];

        auto& p = parameters_[i];
        p.sqrt_gamma_over_mass = std::sqrt(p.gamma * pv.rmass);
        p.accel = pv.force * pv.rmass + this->gen_gaussian_vec(
                  this->noise_coef_ * p.sqrt_gamma_over_mass);
    }
    return;
}

template<typename traitsT>
typename UnderdampedLangevinStepper<traitsT>::real_type
UnderdampedLangevinStepper<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff)
{
    real_type largest_disp2 = real_type(0.0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        auto pv = sys[i]; // particle_view that points i-th particle.
        const auto& param = this->parameters_[i];

        const real_type gamma_dt_over_2           = param.gamma * halfdt_;
        const real_type one_minus_gamma_dt_over_2 = real_type(1.0) - gamma_dt_over_2;

        const auto disp = (dt_ * one_minus_gamma_dt_over_2) * (pv.velocity) +
                halfdt2_ * param.accel;

        pv.position = sys.adjust_position(pv.position + disp);

        pv.velocity *= one_minus_gamma_dt_over_2 *
            (one_minus_gamma_dt_over_2 * one_minus_gamma_dt_over_2 +
             gamma_dt_over_2);
        pv.velocity += (halfdt_ * one_minus_gamma_dt_over_2) * param.accel;

        pv.force = coordinate_type(real_type(0.0), real_type(0.0), real_type(0.0));

        largest_disp2 = std::max(largest_disp2, length_sq(disp));
    }
    sys.largest_displacement() = std::sqrt(largest_disp2);

    // calc f(t+dt)
    ff.calc_force(sys);

    // calc a(t+dt) and v(t+dt), generate noise
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        auto pv = sys[i]; // particle_view that points i-th particle.
        auto& param = this->parameters_[i];
        param.accel = pv.force * pv.rmass + this->gen_gaussian_vec(
                      this->noise_coef_ * param.sqrt_gamma_over_mass);
        pv.velocity += halfdt_ * (real_type(1.0) - param.gamma * halfdt_) * param.accel;
    }

    return time + dt_;
}


} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
