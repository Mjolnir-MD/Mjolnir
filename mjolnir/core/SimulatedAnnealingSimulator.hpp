#ifndef MJOLNIR_CORE_SIMULATED_ANNEALING_SIMULATOR_HPP
#define MJOLNIR_CORE_SIMULATED_ANNEALING_SIMULATOR_HPP
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/ObserverContainer.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>

namespace mjolnir
{

template<typename realT>
struct LinearScheduler
{
    using real_type = realT;

    LinearScheduler(const real_type first, const real_type last)
        : first_(first), last_(last)
    {}
    ~LinearScheduler() = default;

    LinearScheduler(const LinearScheduler&) = default;
    LinearScheduler(LinearScheduler&&)      = default;
    LinearScheduler& operator=(const LinearScheduler&) = default;
    LinearScheduler& operator=(LinearScheduler&&)      = default;

    real_type current(const real_type ratio) const noexcept
    {
        assert(real_type(0.0) <= ratio && ratio <= real_type(1.0));
        return this->first_ * (real_type(1.0) - ratio) + this->last_ * ratio;
    }

  private:
    real_type first_, last_;
};

// XXX : currently, it's not a `CORRECT` simulated annealing because it does not
//       affect to ForceField parameters like DebyeHuckel's debye length.
template<typename traitsT, typename integratorT,
         template<typename> class scheduleT>
class SimulatedAnnealingSimulator final : public SimulatorBase
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using integrator_type = integratorT;
    using system_type     = System<traits_type>;
    using forcefield_type = ForceField<traits_type>;
    using observer_type   = ObserverContainer<traits_type>;
    using scheduler_type  = scheduleT<real_type> ;
    using rng_type        = RandomNumberGenerator<traits_type>;

    SimulatedAnnealingSimulator(
        const std::size_t tstep,     const std::size_t sstep,
        const std::size_t each_step, scheduler_type&&  scheduler,
        system_type&&     sys,       forcefield_type&& ff,
        integrator_type&& integr,    observer_type&&   obs,
        rng_type&& rng)
    : total_step_(tstep), step_count_(0), save_step_(sstep), each_step_(each_step),
      time_(0.0), r_total_step_(1.0 / tstep), scheduler_(scheduler),
      system_(std::move(sys)), ff_(std::move(ff)),
      integrator_(std::move(integr)), observers_(std::move(obs)),
      rng_(std::move(rng))
    {}
    ~SimulatedAnnealingSimulator() override {}

    void initialize() override;
    bool step()       override;
    void run()        override;
    void finalize()   override;

    real_type calc_energy() const {return this->ff_.calc_energy(this->system_);}

    system_type&       system()       noexcept {return system_;}
    system_type const& system() const noexcept {return system_;}

    forcefield_type&       forcefields()       noexcept {return ff_;}
    forcefield_type const& forcefields() const noexcept {return ff_;}

    real_type& time()       noexcept {return time_;}
    real_type  time() const noexcept {return time_;}

  protected:
    std::size_t     total_step_;
    std::size_t     step_count_;
    std::size_t     save_step_;
    std::size_t     each_step_;
    real_type       time_;
    real_type       r_total_step_;
    scheduler_type  scheduler_;
    system_type     system_;
    forcefield_type ff_;
    integrator_type integrator_;
    observer_type   observers_;
    rng_type        rng_;
};

template<typename traitsT, typename integratorT,
         template<typename> class scheduleT>
inline void
SimulatedAnnealingSimulator<traitsT, integratorT, scheduleT>::initialize()
{
    system_.attribute("temperature") = this->scheduler_.current(0);

    this->ff_.initialize(this->system_);
    this->integrator_.initialize(this->system_, this->ff_, this->rng_);

    observers_.initialize(this->total_step_, this->integrator_.delta_t(),
                          this->system_, this->ff_);
    return;
}

template<typename traitsT, typename integratorT,
         template<typename> class scheduleT>
inline bool SimulatedAnnealingSimulator<traitsT, integratorT, scheduleT>::step()
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    if(step_count_ % save_step_ == 0)
    {
        observers_.output(this->step_count_, this->integrator_.delta_t(),
                          this->system_, this->ff_);
    }

    integrator_.step(this->time_, system_, ff_, this->rng_);
    ++step_count_;
    this->time_ = this->step_count_ * integrator_.delta_t();

    if(this->step_count_ % each_step_ == 0)
    {
        MJOLNIR_LOG_SCOPE_DEBUG(if(this->step_count_ % each_step == 0));

        system_.attribute("temperature") =
            this->scheduler_.current(step_count_ * r_total_step_);

        MJOLNIR_LOG_DEBUG("T = ", system_.attribute("temperature"));
    }

    this->integrator_.update(system_);

    return step_count_ < total_step_;
}

template<typename traitsT, typename integratorT, template<typename> class scheduleT>
inline void SimulatedAnnealingSimulator<traitsT, integratorT, scheduleT>::run()
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    for(std::size_t i=0; i<total_step_; ++i)
    {
        this->step();
    }
    assert(this->step_count_ == total_step_);
    return;
}

template<typename traitsT, typename integratorT,
         template<typename> class scheduleT>
inline void
SimulatedAnnealingSimulator<traitsT, integratorT, scheduleT>::finalize()
{
    observers_.output  (this->step_count_, this->integrator_.delta_t(),
                        this->system_, this->ff_);
    observers_.finalize(this->step_count_, this->integrator_.delta_t(),
                        this->system_, this->ff_);
    return;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
namespace mjolnir
{
// BAOAB
extern template class SimulatedAnnealingSimulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >, LinearScheduler>;
extern template class SimulatedAnnealingSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >, LinearScheduler>;
extern template class SimulatedAnnealingSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>, LinearScheduler>;
extern template class SimulatedAnnealingSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>, LinearScheduler>;
// Langevin
extern template class SimulatedAnnealingSimulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >, LinearScheduler>;
extern template class SimulatedAnnealingSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >, LinearScheduler>;
extern template class SimulatedAnnealingSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>, LinearScheduler>;
extern template class SimulatedAnnealingSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>, LinearScheduler>;
}
#endif // SEPARATE_BUILD

#endif // MJOLNIR_SIMULATED_ANNEALING_SIMULATOR_HPP
