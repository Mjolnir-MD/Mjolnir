#ifndef MJOLNIR_OMP_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_OMP_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include <mjolnir/util/aligned_allocator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>

namespace mjolnir
{

// a specialization of UnderdampedLangevinIntegrator for OpenMP implementation.
template<typename realT, template<typename, typename> class boundaryT>
class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = ForceField<traits_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;

  public:

    UnderdampedLangevinIntegrator(const real_type dt,
            std::vector<real_type>&& gamma)
        : dt_(dt), halfdt_(dt / 2), halfdt2_(dt * dt / 2),
          gammas_(std::move(gamma)),
          sqrt_gamma_over_mass_(gammas_.size()),
          acceleration_(gammas_.size())
    {
        assert(this->sqrt_gamma_over_mass_.size() == this->gammas_.size());
        assert(this->acceleration_.size()         == this->gammas_.size());
    }
    ~UnderdampedLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type& rng)
    {
        // initialize temperature and noise intensity
        this->update(sys);

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        }

        ff.calc_force(sys);

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto  rmass = sys.rmass(i);
            const auto& force = sys.force(i);

            sqrt_gamma_over_mass_[i] = std::sqrt(gammas_[i] * rmass);
            acceleration_[i]         = force * rmass +
                this->gen_gaussian_vec(rng, this->noise_coef_ * sqrt_gamma_over_mass_[i]);
        }
        return;
    }

    real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                   rng_type& rng)
    {
        real_type largest_disp2(0.0);

#pragma omp parallel for reduction(max:largest_disp2)
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

            // clear force
            f = math::make_coordinate<coordinate_type>(0, 0, 0);

            largest_disp2 = std::max(largest_disp2, math::length_sq(displacement));
        }

        // This function parallelize itself inside. `parallel` block is not
        // needed here to parallelize it.
        ff.update_margin(2 * std::sqrt(largest_disp2), sys);
        ff.calc_force(sys);

        // calc a(t+dt) and v(t+dt), generate noise
#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto  rm = sys.rmass(i);
            auto&       v  = sys.velocity(i);
            const auto& f  = sys.force(i);
            auto&       a  = this->acceleration_[i];

            a  = f * rm + gen_gaussian_vec(rng, noise_coef_ * sqrt_gamma_over_mass_[i]);
            v += halfdt_ * (1 - gammas_[i] * halfdt_) * a;
        }
        return time + dt_;
    }

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
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
