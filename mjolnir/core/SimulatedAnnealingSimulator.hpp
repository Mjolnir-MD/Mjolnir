#ifndef MJOLNIR_SIMULATED_ANNEALING_SIMULATOR_HPP
#define MJOLNIR_SIMULATED_ANNEALING_SIMULATOR_HPP
#include "SimulatorBase.hpp"
#include "System.hpp"
#include "ForceField.hpp"
#include "Observer.hpp"

namespace mjolnir
{

// TODO: add some options to control temperature
// XXX : currently, it's not a `CORRECT` simulated annealing because it does not
//       affect to ForceField parameters like DebyeHuckel's debye length.
template<typename traitsT, typename integratorT>
class SimulatedAnnealingSimulator final : public SimulatorBase
{
  public:
    typedef traitsT     traits_type;
    typedef integratorT integrator_type;
    typedef System<traits_type>     system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef Observer<traits_type>   observer_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    SimulatedAnnealingSimulator(const std::size_t tstep,
            const real_type   T_first, const real_type T_last,
            system_type&&     sys,     forcefield_type&& ff,
            integrator_type&& integr,  observer_type&&   obs)
    : total_step_(tstep), step_count_(0), time_(0.), r_total_step_(1.0 / tstep),
      first_temperature_(T_first), last_temperature_(T_last),
      system_(std::move(sys)), ff_(std::move(ff)),
      integrator_(std::move(integr)), observer_(std::move(obs))
    {}
    ~SimulatedAnnealingSimulator() override = default;

    void initialize() override;
    bool step()       override;
    void finalize()   override;

    real_type calc_energy() const {return this->ff_.calc_energy(this->system_);}

    system_type&       system()       noexcept {return system_;}
    system_type const& system() const noexcept {return system_;}

    ForceField<traitsT>&       forcefields()       noexcept {return ff_;}
    ForceField<traitsT> const& forcefields() const noexcept {return ff_;}

    real_type& time()       noexcept {return time_;}
    real_type  time() const noexcept {return time_;}

  protected:
    std::size_t     total_step_;
    std::size_t     step_count_;
    real_type       time_;
    real_type       r_total_step_;
    real_type       first_temperature_;
    real_type       last_temperature_;
    system_type     system_;
    forcefield_type ff_;
    integrator_type integrator_;
    observer_type   observer_;
};

template<typename traitsT, typename integratorT>
inline void SimulatedAnnealingSimulator<traitsT, integratorT>::initialize()
{
    this->ff_.initialize(this->system_, integrator_.delta_t());
    this->integrator_.initialize(this->system_, this->ff_);

    observer_.initialize(this->system_, this->ff_);
    observer_.output(0., this->system_, this->ff_);
    return;
}

template<typename traitsT, typename integratorT>
inline bool SimulatedAnnealingSimulator<traitsT, integratorT>::step()
{
    integrator_.step(this->time_, system_, ff_);
    ++step_count_;
    this->time_ = this->step_count_ * integrator_.delta_t();

    // TODO create new object that controls temperature.
    const real_type ratio = step_count_ * r_total_step_;
    system_.attribute("temperature") = this->first_temperature_ * (1 - ratio) +
                                       this->last_temperature_  * ratio;

    this->integrator_.update(system_);

    if(observer_.is_output_time())
    {
        observer_.output(this->time_, this->system_, this->ff_);
    }
    return step_count_ < total_step_;
}

template<typename traitsT, typename integratorT>
inline void SimulatedAnnealingSimulator<traitsT, integratorT>::finalize()
{
    return;
}

} // mjolnir
#endif // MJOLNIR_SIMULATED_ANNEALING_SIMULATOR_HPP
