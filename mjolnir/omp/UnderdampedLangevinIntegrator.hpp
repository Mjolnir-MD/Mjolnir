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
            std::vector<real_type>&& gamma, rng_type&& rng)
        : dt_(dt), halfdt_(dt / 2), halfdt2_(dt * dt / 2),
          rng_(std::move(rng)),
          gammas_(std::move(gamma)),
          sqrt_gamma_over_mass_(gammas_.size()),
          acceleration_(gammas_.size())
    {
        assert(this->sqrt_gamma_over_mass_.size() == this->gammas_.size());
        assert(this->acceleration_.size()         == this->gammas_.size());
    }
    ~UnderdampedLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff)
    {
        // initialize temperature and noise intensity
        this->update(sys);

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        }

#pragma omp parallel
        {
            // calc_force uses `nowait` to speedup. To do that, it needs to be
            // inside a parallel region. So, only for this, we need to wrap
            // `calc_force` with `parallel` region.
            ff.calc_force(sys);
        }

        // merge thread-local forces into master
        sys.merge_forces();

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto  rmass = sys.rmass(i);
            const auto& force = sys.force(i);

            sqrt_gamma_over_mass_[i] = std::sqrt(gammas_[i] * rmass);
            acceleration_[i]         = force * rmass +
                this->gen_gaussian_vec(this->noise_coef_ * sqrt_gamma_over_mass_[i]);
        }
        return;
    }
    real_type step(const real_type time, system_type& sys, forcefield_type& ff)
    {
        real_type largest_disp2(0.0);

#pragma omp parallel for reduction(max:largest_disp2)
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            auto&       p = sys.position(i); // position of i-th particle
            auto&       v = sys.velocity(i); // ditto
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

            largest_disp2 = std::max(largest_disp2, math::length_sq(displacement));
        }

        // This function parallelize itself inside. `parallel` block is not
        // needed here to parallelize it.
        ff.update_margin(2 * std::sqrt(largest_disp2), sys);

#pragma omp parallel
        {
            // calc_force uses `nowait` to speedup. To do that, it needs to be
            // inside a parallel region. So, only for this, we need to wrap
            // `calc_force` with `parallel` region.
            ff.calc_force(sys);
        }

        // merge thread-local forces into master
        sys.merge_forces();

        // calc a(t+dt) and v(t+dt), generate noise
#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto  rm = sys.rmass(i);
            auto&       v  = sys.velocity(i);
            const auto& f  = sys.force(i);
            auto&       a  = this->acceleration_[i];

            a  = f * rm + gen_gaussian_vec(noise_coef_ * sqrt_gamma_over_mass_[i]);
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

    coordinate_type gen_gaussian_vec(const real_type coef)
    {
        return math::make_coordinate<coordinate_type>(
                this->rng_.gaussian(0, coef),
                this->rng_.gaussian(0, coef),
                this->rng_.gaussian(0, coef));
    }

  private:
    real_type dt_;
    real_type halfdt_;
    real_type halfdt2_;
    real_type temperature_;
    real_type noise_coef_;
    rng_type  rng_;

    std::vector<real_type>       gammas_;
    std::vector<real_type>       sqrt_gamma_over_mass_;
    std::vector<coordinate_type> acceleration_;
};

} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
