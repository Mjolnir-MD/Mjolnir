#ifndef MJOLNIR_CORE_BAOAB_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CORE_BAOAB_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// BAOAB Langevin integrator developed by the following papers
// - Leimkuhler B, Matthews C. Appl. Math. Res. Exp. (2013)
// - Leimkuhler B, Matthews C. J. Chem. Phys. (2013)
template<typename traitsT>
class BAOABLangevinIntegrator
{
  public:
    using traits_type     = traitsT;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using matrix33_type   = typename traits_type::matrix33_type;
    using system_type     = System<traitsT>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traitsT>>;
    using rng_type        = RandomNumberGenerator<traits_type>;
    using remover_type    = SystemMotionRemover<traits_type>;
    using variable_key_type = typename system_type::variable_key_type;

  public:

    BAOABLangevinIntegrator(const real_type dt, std::vector<real_type>&& gamma,
                            remover_type&& remover)
        : dt_(dt), halfdt_(dt / 2), gammas_(std::move(gamma)),
          exp_gamma_dt_(gammas_.size()), noise_coeff_ (gammas_.size()),
          remover_(std::move(remover))
    {}
    ~BAOABLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type& rng);

    real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                   rng_type& rng);

    void update(const system_type& sys)
    {
        if(!sys.has_attribute("temperature"))
        {
            throw std::out_of_range("mjolnir::BAOABLangevinIntegrator: "
                "Langevin Integrator requires reference temperature, but "
                "`temperature` is not found in `system.attribute`.");
        }
        this->temperature_ = sys.attribute("temperature");
        this->reset_parameters(sys);
        return;
    }

    real_type delta_t() const noexcept {return dt_;}
    std::vector<real_type> const& parameters() const noexcept {return gammas_;}

  private:

    void reset_parameters(const system_type& sys) noexcept
    {
        const auto kBT = physics::constants<real_type>::kB() * this->temperature_;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto gamma    = this->gammas_.at(i);
            const auto gamma_dt = -1 * gamma * this->dt_;
            this->exp_gamma_dt_.at(i) = std::exp(gamma_dt);
            this->noise_coeff_ .at(i) = std::sqrt(
                    kBT * (1 - std::exp(2 * gamma_dt)) * sys.rmass(i));
        }

        for(const auto& kv : sys.variables())
        {
            const auto& key = kv.first;
            const auto& var = kv.second;

            // force is not initialized yet
            dynvar_params param;
            param.exp_gamma_dt = std::exp(-var.gamma() * this->dt_);
            param.noise_coeff  = std::sqrt(kBT *
                    (real_type(1) - std::exp(-2 * var.gamma() * this->dt_)) /
                    var.m());
            params_for_dynvar_[key] = param;
        }
        return;
    }

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

    remover_type remover_;

    struct dynvar_params
    {
        real_type exp_gamma_dt;
        real_type noise_coeff;
    };
    std::map<variable_key_type, dynvar_params> params_for_dynvar_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own specialization to run it in parallel.
    // So this implementation should not be instanciated with the OpenMP traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT>
void BAOABLangevinIntegrator<traitsT>::initialize(
        system_type& system, forcefield_type& ff, rng_type&)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    if(!ff->constraint().empty())
    {
        MJOLNIR_LOG_WARN("BAOAB langevin integrator does not support constraint"
            " forcefield. [[forcefields.constraint]] will be ignored.");
    }

    // calculate parameters for each particles
    this->update(system);

    // if loaded from MsgPack, we can skip it.
    if( ! system.force_initialized())
    {
        // calculate force
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
typename BAOABLangevinIntegrator<traitsT>::real_type
BAOABLangevinIntegrator<traitsT>::step(
        const real_type time, system_type& sys, forcefield_type& ff, rng_type& rng)
{
    real_type largest_disp2(0.0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto R  = this->gen_R(rng); // random gaussian vector (0 mean, 1 var)
        const auto rm = sys.rmass(i);  // reciprocal mass
        auto&      p  = sys.position(i);
        auto&      v  = sys.velocity(i);
        auto&      f  = sys.force(i);
        const auto expgt = this->exp_gamma_dt_[i]; // exp(- gamma dt)
        coordinate_type dp = math::make_coordinate<coordinate_type>(0, 0, 0);

        v  += this->halfdt_ * rm * f;    // calc v(n+1/3)
        dp += this->halfdt_ * v;         // calc p(n+1/2)
        v  *= expgt;
        v  += this->noise_coeff_[i] * R; // calc v(n+2/3)
        dp += this->halfdt_ * v;         // calc p(n+1)

        // update p(n) -> p(n+1)
        p = sys.adjust_position(p + dp);

        // reset force
        f = math::make_coordinate<coordinate_type>(0, 0, 0);

        // collect largest displacement
        largest_disp2 = std::max(largest_disp2, math::length_sq(dp));
    }
    for(auto& kv : sys.variables()) // assume there are only a few dynvars
    {
        const auto& key = kv.first;
        auto&       var = kv.second;

        const auto& param = params_for_dynvar_.at(key);

        real_type next_v = var.v() + halfdt_ * var.f() / var.m();
        real_type next_x = var.x() + halfdt_ * var.v();
        next_v *= param.exp_gamma_dt;
        next_v += param.noise_coeff * rng.gaussian();
        next_x += halfdt_ * var.v();

        var.update(next_x, next_v, real_type(0));
    }
    sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff->reduce_margin(2 * std::sqrt(largest_disp2), sys);

    // calc f(p(n+1))
    ff->calc_force(sys);

    // calc v(n+2/3) -> v(n+1)
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto rm = sys.rmass(i);  // reciprocal math
        const auto& f = sys.force(i);
        auto&       v = sys.velocity(i);
        v += this->halfdt_ * rm * f;
    }
    for(auto& kv : sys.variables())
    {
        auto& var = kv.second;
        var.update(var.x(), var.v() + this->halfdt_ * var.f() / var.m(), var.f());
    }

    // remove net rotation/translation
    remover_.remove(sys);

    return time + dt_;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_BAOAB_LANGEVIN_INTEGRATOR */
