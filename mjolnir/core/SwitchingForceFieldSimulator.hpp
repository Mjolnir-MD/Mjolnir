#ifndef MJOLNIR_CORE_SWITCHING_FORCEFIELD_SIMULATOR_HPP
#define MJOLNIR_CORE_SWITCHING_FORCEFIELD_SIMULATOR_HPP
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/ObserverContainer.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>

namespace mjolnir
{
// [simulator]
// type = "SwitchingForceField"
// schedule = [
//   # execute from the first one to the last.
//   {until = 10000, forcefield = "open"},  # means [    0, 10000)
//   {until = 20000, forcefield = "close"}, # means [10000, 20000)
//   {until = 30000, forcefield = "open"},  # means [20000, 30000)
//   # step from t = 9999 to 10000, it uses "open".
//   # from t = 10000 to 10001, it uses "close".
// ]
// # ...
// [[forcefields]]
// name = "open"
// [[forcefields.local]]
// # ...
//
// [[forcefields]]
// name = "close"
// [[forcefields.local]]
// # ...
template<typename traitsT, typename integratorT>
class SwitchingForceFieldSimulator final : public SimulatorBase
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using integrator_type = integratorT;
    using system_type     = System<traits_type>;
    using forcefield_type = ForceField<traits_type>;
    using observer_type   = ObserverContainer<traits_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;

  public:

    SwitchingForceFieldSimulator(
        const std::size_t                                  tstep,
        const std::size_t                                  save_step,
        system_type&&                                      sys,
        std::vector<forcefield_type>&&                     ff,
        integrator_type&&                                  integr,
        observer_type&&                                    obs,
        rng_type&&                                         rng,
        std::map<std::string, std::size_t>&&               forcefield_index,
        std::vector<std::pair<std::size_t, std::string>>&& schedule)
        : current_forcefield_(0),
          current_schedule_(0),
          next_switch_step_(0),
          total_step_(tstep),
          step_count_(0),
          save_step_(save_step),
          time_(0.),
          system_(std::move(sys)),
          integrator_(std::move(integr)),
          observers_(std::move(obs)),
          rng_(std::move(rng)),
          forcefields_(std::move(ff)),
          forcefield_index_(std::move(forcefield_index)),
          schedule_(std::move(schedule))
    {}
    ~SwitchingForceFieldSimulator() override {}

    void initialize() override;
    bool step()       override;
    void run()        override;
    void finalize()   override;

    system_type&       system()       noexcept {return system_;}
    system_type const& system() const noexcept {return system_;}

    std::vector<forcefield_type>&       forcefields()       noexcept {return forcefields_;}
    std::vector<forcefield_type> const& forcefields() const noexcept {return forcefields_;}

    real_type& time()       noexcept {return time_;}
    real_type  time() const noexcept {return time_;}

    rng_type&       rng()       noexcept {return rng_;}
    rng_type const& rng() const noexcept {return rng_;}

    // ------------------------------------------------------------------------
    // for testing

    std::map<std::string, std::size_t> const&
    forcefield_index() const noexcept {return forcefield_index_;}

    std::vector<std::pair<std::size_t, std::string>> const&
    schedule() const noexcept {return schedule_;}

  private:

    std::size_t     current_forcefield_; // index of current ff in forcefields_
    std::size_t     current_schedule_;   // current schedule
    std::size_t     next_switch_step_;   // the next step we need to update the above
    std::size_t     total_step_;
    std::size_t     step_count_;
    std::size_t     save_step_;
    real_type       time_;
    system_type     system_;
    integrator_type integrator_;
    observer_type   observers_;
    rng_type        rng_;
    std::vector<forcefield_type>                  forcefields_;
    std::map<std::string, std::size_t>       forcefield_index_;
    std::vector<std::pair<std::size_t, std::string>> schedule_;
};

template<typename traitsT, typename integratorT>
inline void SwitchingForceFieldSimulator<traitsT, integratorT>::initialize()
{
    this->current_schedule_   = 0;

    const auto& sch = schedule_.at(current_schedule_);

    this->next_switch_step_   = sch.first;
    this->current_forcefield_ = this->forcefield_index_.at(sch.second);

    auto& ff = this->forcefields_[this->current_forcefield_];

    this->system_.initialize(this->rng_);
    ff.initialize(this->system_);
    this->integrator_.initialize(this->system_, ff, this->rng_);

    observers_.initialize(this->total_step_, this->integrator_.delta_t(),
                          this->system_, ff);

    return;
}

template<typename traitsT, typename integratorT>
inline bool SwitchingForceFieldSimulator<traitsT, integratorT>::step()
{
    if(step_count_ % save_step_ == 0)
    {
        observers_.output(this->step_count_, this->integrator_.delta_t(),
                          this->system_, forcefields_[current_forcefield_]);
    }

    integrator_.step(this->time_, system_, forcefields_[current_forcefield_],
                     this->rng_);

    ++step_count_;
    this->time_ = this->step_count_ * integrator_.delta_t();

    // the step reaches to the `until`. change forcefield.
    // if until == total_step, the next schedule does not exists.
    // so first check the simulation continues.
    if(step_count_ < total_step_ && step_count_ == next_switch_step_)
    {
        this->current_schedule_  += 1;
        const auto& sch = schedule_.at(current_schedule_);

        this->next_switch_step_   = sch.first;
        this->current_forcefield_ = forcefield_index_.at(sch.second);

        // initialize spatial partition (e.g. cell lists) and system.topology
        forcefields_[current_forcefield_].initialize(this->system_);
        // initialize forces with the current forcefield. the previous forces
        // will be zero-cleared. but velocities are kept.
        integrator_.initialize(this->system_, forcefields_[current_forcefield_],
                               this->rng_);
        // Observers need to be updated because forcefield changed.
        // especially, EnergyObserver needed to be updated.
        observers_.update(this->step_count_, this->integrator_.delta_t(),
                          this->system_, forcefields_[current_forcefield_]);
    }
    return step_count_ < total_step_;
}

template<typename traitsT, typename integratorT>
inline void SwitchingForceFieldSimulator<traitsT, integratorT>::run()
{
    for(std::size_t i=0; i<total_step_; ++i)
    {
        this->step();
    }
    assert(this->step_count_ == total_step_);
    return;
}

template<typename traitsT, typename integratorT>
inline void SwitchingForceFieldSimulator<traitsT, integratorT>::finalize()
{
    observers_.output(this->step_count_, this->integrator_.delta_t(),
                      this->system_, forcefields_[current_forcefield_]);
    observers_.finalize(this->total_step_, this->integrator_.delta_t(),
                        this->system_, forcefields_[current_forcefield_]);
    return;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BAOABLangevinIntegrator.hpp>
#include <mjolnir/core/UnderdampedLangevinIntegrator.hpp>
#include <mjolnir/core/VelocityVerletIntegrator.hpp>
namespace mjolnir
{
// BAOAB
extern template class SwitchingForceFieldSimulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// Langevin
extern template class SwitchingForceFieldSimulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// VelVerlet
extern template class SwitchingForceFieldSimulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class SwitchingForceFieldSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
}
#endif // SEPARATE_BUILD


#endif// MJOLNIR_CORE_SWITCHING_FORCEFIELD_SIMULATOR_HPP
