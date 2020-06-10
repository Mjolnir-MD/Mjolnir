#ifndef MJOLNIR_CORE_ENERGY_CALCULATION_SIMULATOR_HPP
#define MJOLNIR_CORE_ENERGY_CALCULATION_SIMULATOR_HPP
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/ObserverContainer.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/LoaderBase.hpp>
#include <mjolnir/core/XYZLoader.hpp>
#include <mjolnir/core/DCDLoader.hpp>
#include <mjolnir/core/TRRLoader.hpp>

namespace mjolnir
{

template<typename traitsT>
class EnergyCalculationSimulator final : public SimulatorBase
{
  public:
    using traits_type      = traitsT;
    using real_type        = typename traits_type::real_type;
    using coordinate_type  = typename traits_type::coordinate_type;
    using system_type      = System<traits_type>;
    using forcefield_type  = std::unique_ptr<ForceFieldBase<traits_type>>;
    using observer_type    = ObserverContainer<traits_type>;
    using loader_base_type = LoaderBase<traits_type>;
    using loader_type      = std::unique_ptr<loader_base_type>;

  public:

    EnergyCalculationSimulator(const std::size_t total_step, loader_type&& ld,
            system_type&& sys, forcefield_type&& ff, observer_type&& obs)
        : step_count_(0),     total_step_(total_step),
          ld_(std::move(ld)), sys_(std::move(sys)),
          ff_(std::move(ff)), obs_(std::move(obs))
    {}
    ~EnergyCalculationSimulator() override {}

    void initialize() override;
    bool step()       override;
    void run()        override;
    void finalize()   override;

    system_type&       system()       noexcept {return sys_;}
    system_type const& system() const noexcept {return sys_;}

    forcefield_type&       forcefields()       noexcept {return ff_;}
    forcefield_type const& forcefields() const noexcept {return ff_;}

    loader_type&       loader()       noexcept {return ld_;}
    loader_type const& loader() const noexcept {return ld_;}

  protected:

    std::size_t     step_count_, total_step_;
    loader_type     ld_;
    system_type     sys_;
    forcefield_type ff_;
    observer_type   obs_;
};

template<typename traitsT>
inline void EnergyCalculationSimulator<traitsT>::initialize()
{
    // instead of initializing simulator (it generate velocities), read the
    // current state from the loader.
    //
    // To obtain the total number of snapshots, the loader is already initialized.
    ld_->load_next(this->sys_);

    // sys.initialize() generates the initial velocity, so here it is not needed

    this->ff_->initialize(this->sys_);

    this->obs_.initialize(this->total_step_, 0.0, this->sys_, this->ff_);
    return;
}

template<typename traitsT>
inline bool EnergyCalculationSimulator<traitsT>::step()
{
    // this calculates the energy.
    obs_.output(this->step_count_, 0.0, this->sys_, this->ff_);

    if(!ld_->load_next(sys_))
    {
        // if loader returns false, it means it fails to read the next traj.
        return false;
    }

    ff_->update(sys_); // force update the neighboring lists

    ++step_count_;

    return step_count_ < total_step_;
}

template<typename traitsT>
inline void EnergyCalculationSimulator<traitsT>::run()
{
    while(this->step()) {};
    return;
}

template<typename traitsT>
inline void EnergyCalculationSimulator<traitsT>::finalize()
{
    obs_.finalize(this->total_step_, 0.0, this->sys_, this->ff_);
    return;
}

} // mjolnir
#endif// MJOLNIR_CORE_ENERGY_CALCULATION_SIMULATOR_HPP
