#ifndef MJOLNIR_CORE_MOLECULAR_DYNAMICS_SIMULATOR_HPP
#define MJOLNIR_CORE_MOLECULAR_DYNAMICS_SIMULATOR_HPP
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/ObserverContainer.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/core/MsgPackSaver.hpp>

namespace mjolnir
{

template<typename traitsT, typename integratorT>
class MolecularDynamicsSimulator final : public SimulatorBase
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using integrator_type = integratorT;
    using system_type     = System<traits_type>;
    using forcefield_type = std::unique_ptr<ForceFieldBase<traits_type>>;
    using observer_type   = ObserverContainer<traits_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;
    using saver_type      = MsgPackSaver<traits_type>;

    MolecularDynamicsSimulator(const std::size_t tstep,
        const std::size_t save_step, const std::size_t checkpoint,
        system_type&& sys, forcefield_type&& ff,
        integrator_type&& integr, observer_type&& obs,
        rng_type&& rng)
        : total_step_(tstep), step_count_(0), save_step_(save_step),
          checkpoint_(checkpoint), time_(0),
          system_(std::move(sys)), ff_(std::move(ff)),
          integrator_(std::move(integr)), observers_(std::move(obs)),
          saver_(observers_.prefix()), rng_(std::move(rng))
    {}
    ~MolecularDynamicsSimulator() override {}

    void initialize() override;
    bool step()       override;
    void run()        override;
    void finalize()   override;

    system_type&       system()       noexcept {return system_;}
    system_type const& system() const noexcept {return system_;}

    forcefield_type&       forcefields()       noexcept {return ff_;}
    forcefield_type const& forcefields() const noexcept {return ff_;}

    real_type& time()       noexcept {return time_;}
    real_type  time() const noexcept {return time_;}

    rng_type&       rng()       noexcept {return rng_;}
    rng_type const& rng() const noexcept {return rng_;}

  protected:
    std::size_t     total_step_;
    std::size_t     step_count_;
    std::size_t     save_step_;
    std::size_t     checkpoint_;
    real_type       time_;
    system_type     system_;
    forcefield_type ff_;
    integrator_type integrator_;
    observer_type   observers_;
    saver_type      saver_;
    rng_type        rng_;
};

template<typename traitsT, typename integratorT>
void MolecularDynamicsSimulator<traitsT, integratorT>::initialize()
{
    this->system_.initialize(this->rng_);
    this->ff_->initialize(this->system_);
    this->integrator_.initialize(this->system_, this->ff_, this->rng_);

    observers_.initialize(this->total_step_, this->save_step_,
                          this->integrator_.delta_t(), this->system_, this->ff_);
    return;
}

template<typename traitsT, typename integratorT>
bool MolecularDynamicsSimulator<traitsT, integratorT>::step()
{
    if(step_count_ % save_step_ == 0)
    {
        observers_.output(this->step_count_, this->integrator_.delta_t(),
                          this->system_, this->ff_);
    }
    if(step_count_ % checkpoint_ == 0)
    {
        saver_.save(this->system_);
        saver_.save(this->rng_);
    }

    integrator_.step(this->time_, system_, ff_, this->rng_);
    ++step_count_;
    this->time_ = this->step_count_ * integrator_.delta_t();

    return step_count_ < total_step_;
}

template<typename traitsT, typename integratorT>
void MolecularDynamicsSimulator<traitsT, integratorT>::run()
{
    for(std::size_t i=0; i<total_step_; ++i)
    {
        this->step();
    }
    assert(this->step_count_ == total_step_);
    return;
}

template<typename traitsT, typename integratorT>
void MolecularDynamicsSimulator<traitsT, integratorT>::finalize()
{
    observers_.output(this->step_count_, this->integrator_.delta_t(),
                      this->system_, this->ff_);
    observers_.finalize(this->total_step_, this->integrator_.delta_t(),
                        this->system_, this->ff_);
    saver_.save(this->system_);
    saver_.save(this->rng_);
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
extern template class MolecularDynamicsSimulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// Langevin
extern template class MolecularDynamicsSimulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
// VelVerlet
extern template class MolecularDynamicsSimulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>;
extern template class MolecularDynamicsSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>;
}
#endif // SEPARATE_BUILD

#endif /* MJOLNIR_MOLECULAR_DYNAMICS_SIMULATOR */
