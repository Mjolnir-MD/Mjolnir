#ifndef MJOLNIR_MOLECULAR_DYNAMICS_SIMULATOR
#define MJOLNIR_MOLECULAR_DYNAMICS_SIMULATOR
#include "SimulatorBase.hpp"
#include "System.hpp"
#include "ForceField.hpp"
#include "Observer.hpp"

namespace mjolnir
{

template<typename traitsT, typename integratorT>
class MolecularDynamicsSimulator final : public SimulatorBase
{
  public:
    typedef traitsT     traits_type;
    typedef integratorT integrator_type;
    typedef System<traits_type>     system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef Observer<traits_type>   observer_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    MolecularDynamicsSimulator(const std::size_t tstep, const std::size_t save_step,
                system_type&& sys, forcefield_type&& ff,
                integrator_type&& integr, observer_type&& obs)
    : total_step_(tstep), step_count_(0), save_step_(save_step), time_(0.),
      system_(std::move(sys)), ff_(std::move(ff)),
      integrator_(std::move(integr)), observer_(std::move(obs))
    {}
    ~MolecularDynamicsSimulator() override = default;

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
    real_type       time_;
    system_type     system_;
    forcefield_type ff_;
    integrator_type integrator_;
    observer_type   observer_;
};

template<typename traitsT, typename integratorT>
inline void MolecularDynamicsSimulator<traitsT, integratorT>::initialize()
{
    this->ff_.initialize(this->system_);
    this->integrator_.initialize(this->system_, this->ff_);

    observer_.initialize(this->system_, this->ff_);
    observer_.output(0., this->system_, this->ff_);
    return;
}

template<typename traitsT, typename integratorT>
inline bool MolecularDynamicsSimulator<traitsT, integratorT>::step()
{
    integrator_.step(this->time_, system_, ff_);
    ++step_count_;
    this->time_ = this->step_count_ * integrator_.delta_t();

    if(step_count_ % save_step_ == 0)
    {
        observer_.output(this->time_, this->system_, this->ff_);
    }
    return step_count_ < total_step_;
}

template<typename traitsT, typename integratorT>
inline void MolecularDynamicsSimulator<traitsT, integratorT>::finalize()
{
//     observer_.output(time_, this->system_);
    return;
}

} // mjolnir
#endif /* MJOLNIR_MOLECULAR_DYNAMICS_SIMULATOR */
