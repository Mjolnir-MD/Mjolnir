#ifndef MJOLNIR_CORE_GJ_F_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_CORE_GJ_F_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/SystemMotionRemover.hpp>
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

// "A simple and effective Verlet-type algorithm for simulating Langevin dynamics"
// Niels Gronbech-Jensen and Oded Farago, Molecular Physics, 2013, Vol.111 No.8, 983-991

template<typename traitsT>
class GJFNVTLangevinIntegrator
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

    GJFNVTLangevinIntegrator(const real_type dt, std::vector<real_type>&& alpha,
                            remover_type&& remover)
        : dt_(dt), halfdt_(dt / 2), alphas_(std::move(alpha)),
          betas_(alphas_.size()),
          bs_(alphas_.size()),
          vel_coefs_(alphas_.size()),
          noise_(alphas_.size()),
          remover_(std::move(remover))
    {}
    ~GJFNVTLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type&)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if( ! ff->constraint().empty())
        {
            MJOLNIR_LOG_WARN("G-JF langevin integrator does not support"
                "  constraints. [[forcefields.constraint]] will be ignored.");
        }

        // calculate parameters for each particles
        this->update(sys);

        // if loaded from MsgPack, we can skip it.
        if( ! sys.force_initialized())
        {
            // calculate force
            for(std::size_t i=0; i<sys.size(); ++i)
            {
                sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
            }
            sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);
            ff->calc_force(sys);
        }
        return;
    }

    real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                   rng_type& rng);

    void update(const system_type& sys)
    {
        if(!sys.has_attribute("temperature"))
        {
            throw std::out_of_range("mjolnir::GJFNVTLangevinIntegrator: "
                "Langevin Integrator requires reference temperature, but "
                "`temperature` is not found in `system.attribute`.");
        }
        this->temperature_ = sys.attribute("temperature");
        this->reset_parameters(sys);
        return;
    }

    real_type delta_t() const noexcept {return dt_;}
    std::vector<real_type> const& parameters() const noexcept {return alphas_;}

  private:

    void reset_parameters(const system_type& sys) noexcept
    {
        alphas_   .resize(sys.size());
        betas_    .resize(sys.size());
        bs_       .resize(sys.size());
        vel_coefs_.resize(sys.size());
        noise_    .resize(sys.size());
        const auto kBT = physics::constants<real_type>::kB() * this->temperature_;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto alpha   = this->alphas_.at(i);
            const auto m       = sys.mass(i);
            this->bs_.at(i)    = 1.0 / (1.0 + alpha * dt_ * 0.5 / m);
            this->betas_.at(i) = std::sqrt(2 * alpha * kBT * dt_);
            this->vel_coefs_.at(i) = 1.0 - alpha * bs_.at(i) * dt_ / m;
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

    std::vector<real_type> alphas_;
    std::vector<real_type> betas_;
    std::vector<real_type> bs_;
    std::vector<real_type> vel_coefs_;
    std::vector<coordinate_type> noise_;

    remover_type remover_;
};

template<typename traitsT>
typename GJFNVTLangevinIntegrator<traitsT>::real_type
GJFNVTLangevinIntegrator<traitsT>::step(const real_type time,
        system_type& sys, forcefield_type& ff, rng_type& rng)
{
    real_type largest_disp2(0);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto rm = sys.rmass(i);  // reciprocal mass
        auto&       p = sys.position(i);
        auto&       v = sys.velocity(i);
        auto&       f = sys.force(i);
        auto&    beta = this->noise_[i];
        const auto& b = this->bs_[i];

        // beta^(n+1) = N(0, sqrt(2alpha kBT dt))
        // v^(n+1/2)  = v^(n) + dt/2m f^(n) + beta^(n+1) / 2m
        // r^(n+1)    = r^(n) + b dt v^(n+1/2)
        // f^(n+1)    = -dU/dr|r^(n+1)
        // v^(n+1)    = (1 - alpha*b*dt/m) * v^(n+1/2) + dt/2m f^(n+1) + beta^(n+1) / 2m

        beta = this->gen_R(rng) * betas_[i] * rm;
        v   += ((dt_ * rm) * f + beta) * real_type(0.5);

        const auto dp = (b * dt_) * v;

        p   += dp;
        p    = sys.adjust_position(p);
        v   *= vel_coefs_[i];

        // reset force
        f = math::make_coordinate<coordinate_type>(0, 0, 0);

        // collect largest displacement
        largest_disp2 = std::max(largest_disp2, math::length_sq(dp));
    }
    sys.virial() = matrix33_type(0,0,0, 0,0,0, 0,0,0);

    // update neighbor list; reduce margin, reconstruct the list if needed
    ff->reduce_margin(2 * std::sqrt(largest_disp2), sys);
    ff->calc_force(sys);

    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto  rm   = sys.rmass(i);  // reciprocal math
        const auto& f    = sys.force(i);
        auto&       v    = sys.velocity(i);
        const auto& beta = this->noise_[i];

        v += ((dt_ * rm) * f + beta) * real_type(0.5);
    }

    // remove net rotation/translation
    remover_.remove(sys);

    return time + dt_;
}

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GJFNVTLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class GJFNVTLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class GJFNVTLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GJFNVTLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif// MJOLNIR_CORE_GJ_F_LANGEVIN_INTEGRATOR_HPP
