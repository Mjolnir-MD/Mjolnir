#ifndef MJOLNIR_OMP_BAOAB_LANGEVIN_INTEGRATOR_HPP
#define MJOLNIR_OMP_BAOAB_LANGEVIN_INTEGRATOR_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/omp/ForceField.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/omp/SystemMotionRemover.hpp>
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>

namespace mjolnir
{

// a specialization of BAOAB Langevin integrator for OpenMP implementation
template<typename realT, template<typename, typename> class boundaryT>
class BAOABLangevinIntegrator<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using boundary_type   = typename traits_type::boundary_type;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = ForceField<traits_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;
    using remover_type    = SystemMotionRemover<traits_type>;

  public:

    BAOABLangevinIntegrator(const real_type dt,
            std::vector<real_type>&& gamma, remover_type&& remover)
        : dt_(dt), halfdt_(dt / 2),
          gammas_(std::move(gamma)),
          exp_gamma_dt_(gammas_.size()),
          noise_coeff_ (gammas_.size()),
          remover_(std::move(remover))
    {}
    ~BAOABLangevinIntegrator() = default;

    void initialize(system_type& sys, forcefield_type& ff, rng_type&)
    {
        // calculate parameters for each particles
        this->update(sys);

#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        ff.calc_force(sys);
        return;
    }

    real_type step(const real_type time, system_type& sys, forcefield_type& ff,
                   rng_type& rng)
    {
        real_type largest_disp2(0.0);

#pragma omp parallel for reduction(max:largest_disp2)
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

        // ------------------------------------------------------------------
        // update neighbor list; reduce margin, reconstruct the list if needed
        ff.reduce_margin(2 * std::sqrt(largest_disp2), sys);

        // calc f(p(n+1))
        ff.calc_force(sys);

        // ------------------------------------------------------------------
        // calc v(n+2/3) -> v(n+1)
#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto rm = sys.rmass(i);  // reciprocal math
            const auto& f = sys.force(i);
            auto&       v = sys.velocity(i);
            v += this->halfdt_ * rm * f;
        }

        remover_.remove(sys);

        return time + dt_;
    }

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
#pragma omp parallel for
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            const auto gamma    = this->gammas_.at(i);
            const auto gamma_dt = -1 * gamma * this->dt_;
            this->exp_gamma_dt_.at(i) = std::exp(gamma_dt);
            this->noise_coeff_ .at(i) = std::sqrt(
                    kBT * (1 - std::exp(2 * gamma_dt)) * sys.rmass(i));
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
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
extern template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
extern template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class BAOABLangevinIntegrator<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_BAOAB_LANGEVIN_INTEGRATOR */
