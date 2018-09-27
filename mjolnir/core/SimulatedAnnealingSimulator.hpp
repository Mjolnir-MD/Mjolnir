#ifndef MJOLNIR_SIMULATED_ANNEALING_SIMULATOR_HPP
#define MJOLNIR_SIMULATED_ANNEALING_SIMULATOR_HPP
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/Observer.hpp>

namespace mjolnir
{

template<typename realT>
struct linear_schedule
{
    using real_type = realT;

    linear_schedule(const real_type first, const real_type last)
        : first_(first), last_(last)
    {}
    ~linear_schedule() = default;

    linear_schedule(const linear_schedule&) = default;
    linear_schedule(linear_schedule&&)      = default;
    linear_schedule& operator=(const linear_schedule&) = default;
    linear_schedule& operator=(linear_schedule&&)      = default;

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
    typedef traitsT     traits_type;
    typedef integratorT integrator_type;
    typedef System<traits_type>     system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef Observer<traits_type>   observer_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef scheduleT<real_type> scheduler_type;
    typedef progress_bar<50> progress_bar_type;

    SimulatedAnnealingSimulator(
        const std::size_t tstep,     const std::size_t sstep,
        const std::size_t each_step, scheduler_type&&  scheduler,
        system_type&&     sys,       forcefield_type&& ff,
        integrator_type&& integr,    observer_type&&   obs)
    : total_step_(tstep), step_count_(0), save_step_(sstep),
      time_(0.0), r_total_step_(1.0 / tstep),
      each_step_(each_step), scheduler_(scheduler),
      system_(std::move(sys)), ff_(std::move(ff)),
      integrator_(std::move(integr)), observer_(std::move(obs)),
      progress_bar_(tstep)
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
    std::size_t     save_step_;
    std::size_t     each_step_;
    real_type       time_;
    real_type       r_total_step_;
    scheduler_type  scheduler_;
    system_type     system_;
    forcefield_type ff_;
    integrator_type integrator_;
    observer_type   observer_;
    progress_bar_type progress_bar_;
};

template<typename traitsT, typename integratorT,
         template<typename> class scheduleT>
inline void
SimulatedAnnealingSimulator<traitsT, integratorT, scheduleT>::initialize()
{
    system_.attribute("temperature") = this->scheduler_.current(0);

    this->ff_.initialize(this->system_);
    this->integrator_.initialize(this->system_, this->ff_);

    observer_.initialize(this->system_, this->ff_);
    observer_.output(0., this->system_, this->ff_);
    return;
}

template<typename traitsT, typename integratorT,
         template<typename> class scheduleT>
inline bool SimulatedAnnealingSimulator<traitsT, integratorT, scheduleT>::step()
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_SCOPE_DEBUG(SimulatedAnnealingSimulator::step, 0);

    integrator_.step(this->time_, system_, ff_);
    ++step_count_;
    this->time_ = this->step_count_ * integrator_.delta_t();

    if(this->step_count_ % each_step_ == 0)
    {
        MJOLNIR_SCOPE_DEBUG(if(this->step_count_ % each_step == 0), 1);

        system_.attribute("temperature") =
            this->scheduler_.current(step_count_ * r_total_step_);

        MJOLNIR_LOG_DEBUG("T = ", system_.attribute("temperature"));
    }

    this->integrator_.update(system_);

    if(step_count_ % save_step_ == 0)
    {
        observer_.output(this->time_, this->system_, this->ff_);
        std::cerr << progress_bar_.format(this->step_count_);
    }
    return step_count_ < total_step_;
}

template<typename traitsT, typename integratorT,
         template<typename> class scheduleT>
inline void
SimulatedAnnealingSimulator<traitsT, integratorT, scheduleT>::finalize()
{
    std::cerr << progress_bar_.format(this->total_step_) << std::endl;
    return;
}

} // mjolnir
#endif // MJOLNIR_SIMULATED_ANNEALING_SIMULATOR_HPP
