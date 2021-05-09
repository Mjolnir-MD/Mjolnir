#ifndef MJOLNIR_CORE_UNDERDAMPED_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CORE_UNDERDAMPED_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// Simple underdamped Langevin integrator.
// The implementation is based on the article written by J. D. Honeycutt and
// D. Thirumalai (1992) Biopolymers and also the article written by Z. Guo and
// D. Thirumalai (1995) Biopolymers. Same algorithm used as default Langevin
// integrator in CafeMol that is developed by H. Kenzaki et al., (2011) JCTC.
template<typename traitsT>
class UnderdampedLangevinIntegrator
{
  public:
    using traits_type     = traitsT;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traitsT>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traitsT>>;
    using rng_type        = RandomNumberGenerator<traits_type>;
    using remover_type    = SystemMotionRemover<traits_type>;

  public:

    UnderdampedLangevinIntegrator(
            const real_type dt, std::vector<real_type>&& gamma,
            remover_type&& remover)
        : dt_(dt), halfdt_(dt / 2), halfdt2_(dt * dt / 2),
          sqrt_gamma_over_mass_(gamma.size()),
          acceleration_(gamma.size()),
          remover_(std::move(remover))
    {
        gammas_ = std::move(gamma);
    }
    ~UnderdampedLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type& rng);
    real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                   rng_type& rng);

    real_type delta_t() const noexcept {return dt_;}
    void  set_delta_t(const real_type dt) noexcept
    {
        dt_ = dt; halfdt_ = dt / 2; halfdt2_ = dt * dt / 2;
    }

    void update(const system_type& sys)
    {
        if(!sys.has_attribute("temperature"))
        {
            throw std::out_of_range("mjolnir::UnderdampedLangevinIntegrator: "
                "Langevin Integrator requires reference temperature, but "
                "`temperature` is not found in `system.attribute`.");
        }

        this->temperature_ = sys.attribute("temperature");
        this->noise_coef_  = std::sqrt(2 * physics::constants<real_type>::kB() *
                                       this->temperature_ / dt_);
    }

    std::vector<real_type> const& parameters() const noexcept {return gammas_;}

  private:

    coordinate_type gen_gaussian_vec(rng_type& rng, const real_type coef)
    {
        return math::make_coordinate<coordinate_type>(
                rng.gaussian(0, coef),
                rng.gaussian(0, coef),
                rng.gaussian(0, coef));
    }

  private:
    real_type dt_;
    real_type halfdt_;
    real_type halfdt2_;
    real_type temperature_;
    real_type noise_coef_;

    std::vector<real_type>       gammas_;
    std::vector<real_type>       sqrt_gamma_over_mass_;
    std::vector<coordinate_type> acceleration_;

    remover_type remover_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own specialization to run it in parallel.
    // So this implementation should not be instanciated with the OpenMP traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT>
void UnderdampedLangevinIntegrator<traitsT>::initialize(
        system_type& system, forcefield_type& ff, rng_type& rng)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if( ! ff->constraint().empty())
    {
        MJOLNIR_LOG_WARN("Underdamped langevin integrator does not support "
            "constraint forcefield. [[forcefields.constraint]] will be ignored.");
    }

    // initialize temperature and noise intensity
    this->update(system);

    // if loaded from MsgPack, we can skip it.
    if( ! system.force_initialized())
    {
        for(std::size_t i=0; i<system.size(); ++i)
        {
            system.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

        // calculate force
        ff->calc_force(system);

        for(std::size_t i=0; i<system.size(); ++i)
        {
            const auto  rmass = system.rmass(i);
            const auto& force = system.force(i);

            sqrt_gamma_over_mass_[i] = std::sqrt(gammas_[i] * rmass);
            acceleration_[i] = force * rmass +
                this->gen_gaussian_vec(rng, this->noise_coef_ * sqrt_gamma_over_mass_[i]);
        }
    }
    return;
}

template<typename traitsT>
typename UnderdampedLangevinIntegrator<traitsT>::real_type
UnderdampedLangevinIntegrator<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    real_type largest_disp2(0.0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        auto&       p = sys.position(i); // position of i-th particle
        auto&       v = sys.velocity(i); // ditto
        auto&       f = sys.force(i);    // ditto
        const auto& a = this->acceleration_[i];
        auto    gamma = this->gammas_[i];

        const real_type gamma_dt_over_2           = gamma * halfdt_;
        const real_type one_minus_gamma_dt_over_2 = 1 - gamma_dt_over_2;

        const auto displacement =
            (dt_ * one_minus_gamma_dt_over_2) * v + halfdt2_ * a;

        p = sys.adjust_position(p + displacement);

        v *= one_minus_gamma_dt_over_2 * (one_minus_gamma_dt_over_2 *
             one_minus_gamma_dt_over_2 + gamma_dt_over_2);
        v += (halfdt_ * one_minus_gamma_dt_over_2) * a;

        f = math::make_coordinate<coordinate_type>(0, 0, 0);

        largest_disp2 = std::max(largest_disp2, math::length_sq(displacement));
    }
    sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff->reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(t+dt)
    ff->calc_force(sys);

    // calc a(t+dt) and v(t+dt), generate noise
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto  rm = sys.rmass(i);
        auto&       v  = sys.velocity(i);
        const auto& f  = sys.force(i);
        auto&       a  = this->acceleration_[i];

        a  = f * rm + gen_gaussian_vec(rng, noise_coef_ * sqrt_gamma_over_mass_[i]);
        v += halfdt_ * (1 - gammas_[i] * halfdt_) * a;
    }

    // remove net rotation/translation
    remover_.remove(sys);

    return time + dt_;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
