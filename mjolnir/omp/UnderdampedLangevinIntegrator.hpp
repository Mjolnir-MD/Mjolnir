#ifndef MJOLNIR_OMP_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_OMP_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/SystemMotionRemover.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>

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
    using forcefield_type = std::unique_ptr<ForceFieldBase<traits_type>>;
    using rng_type        = RandomNumberGenerator<traits_type>;
    using remover_type    = SystemMotionRemover<traits_type>;
    using variable_key_type = typename system_type::variable_key_type;

  public:

    UnderdampedLangevinIntegrator(const real_type dt,
            std::vector<real_type>&& gamma, remover_type&& remover)
        : dt_(dt), halfdt_(dt / 2), halfdt2_(dt * dt / 2),
          gammas_(std::move(gamma)),
          sqrt_gamma_over_mass_(gammas_.size()),
          acceleration_(gammas_.size()),
          remover_(std::move(remover))
    {
        assert(this->sqrt_gamma_over_mass_.size() == this->gammas_.size());
        assert(this->acceleration_.size()         == this->gammas_.size());
    }
    ~UnderdampedLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type& rng)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        if(!ff->constraint().empty())
        {
            MJOLNIR_LOG_WARN("UnderdampedLangevin integrator does not support "
                "constraints. [[forcefields.constraint]] will be ignored.");
        }

        // initialize temperature and noise intensity
        this->update(sys);

        // if loaded from MsgPack, we can skip it.
        if( ! sys.force_initialized())
        {
#pragma omp parallel for
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
            }
            for(auto& kv : sys.variables())
            {
                auto& var = kv.second;
                var.update(var.x(), var.v(), real_type(0));
            }
            ff->calc_force(sys);
        }

        // here we can assume that forces are initialized.

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto  rmass = sys.rmass(i);
            const auto& force = sys.force(i);

            acceleration_[i] = force * rmass +
                this->gen_gaussian_vec(rng, this->noise_coef_ * sqrt_gamma_over_mass_[i]);
        }
        for(const auto& kv : sys.variables())
        {
            const auto& key = kv.first;
            const auto& var = kv.second;
            auto& param = params_for_dynvar_.at(key);

            param.accel = var.f() / var.m() +
                this->noise_coef_ * param.sqrt_gamma_over_mass * rng.gaussian();
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
        for(auto& kv : sys.variables())
        {
            const auto& key = kv.first;
            auto&       var = kv.second;

            const auto& param = params_for_dynvar_.at(key);

            const real_type one_minus_gamma_dt_over_2 = real_type(1) - param.gamma_dt_over_2;

            const auto next_x = var.x() + dt_ * one_minus_gamma_dt_over_2 * var.v() + halfdt2_ * param.accel;
            const auto next_v = var.v() * one_minus_gamma_dt_over_2 * (one_minus_gamma_dt_over_2 *
                                          one_minus_gamma_dt_over_2 + param.gamma_dt_over_2) +
                               (halfdt_ * one_minus_gamma_dt_over_2) * param.accel;

            var.update(next_x, next_v, real_type(0));
        }

        // This function parallelize itself inside. `parallel` block is not
        // needed here to parallelize it.
        ff->reduce_margin(2 * std::sqrt(largest_disp2), sys);
        ff->calc_force(sys);

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

        for(auto& kv : sys.variables())
        {
            const auto& key = kv.first;
            auto&       var = kv.second;

            auto& param = params_for_dynvar_.at(key);
            param.accel = var.f() / var.m() +
                noise_coef_ * param.sqrt_gamma_over_mass * rng.gaussian();

            const auto next_v = var.v() + halfdt_ * (1 - var.gamma() * halfdt_) * param.accel;
            var.update(var.x(), next_v, var.f());
        }

        remover_.remove(sys);

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

        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sqrt_gamma_over_mass_[i] = std::sqrt(gammas_[i] * sys.rmass(i));
        }

        for(const auto& kv : sys.variables())
        {
            const auto& key = kv.first;
            const auto& var = kv.second;

            // force is not initialized yet
            dynvar_params param;
            param.accel                = real_type(0);
            param.gamma_dt_over_2      = var.gamma() * halfdt_;
            param.sqrt_gamma_over_mass = std::sqrt(var.gamma() / var.m());
            params_for_dynvar_[key] = param;
        }
    }
    void update(const system_type& sys, const real_type dt) noexcept
    {
        this->dt_ = dt;
        this->halfdt_ = dt / 2;
        this->halfdt2_ = dt * dt / 2;
        this->update(sys);
        return;
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

    // {acceleration, gamma dt / 2, 1 - gamma dt / 2, sqrt(gamma / m)}
    struct dynvar_params
    {
        real_type accel;
        real_type gamma_dt_over_2;
        real_type sqrt_gamma_over_mass;
    };
    std::map<variable_key_type, dynvar_params> params_for_dynvar_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class UnderdampedLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
