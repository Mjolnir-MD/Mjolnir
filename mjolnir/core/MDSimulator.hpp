#ifndef MJOLNIR_MD_SIMULATOR
#define MJOLNIR_MD_SIMULATOR
#include "SimulatorBase.hpp"
#include "System.hpp"
#include "ForceField.hpp"
#include "Observer.hpp"

namespace mjolnir
{

template<typename traitsT, typename integratorT>
class MDSimulator final : public SimulatorBase
{
  public:
    typedef traitsT     traits_type;
    typedef integratorT integrator_type;
    typedef System<traits_type>     system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef Observer<traits_type>   observer_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    MDSimulator(const std::size_t tstep, system_type&& sys, forcefield_type&& ff,
                integrator_type&& integr, observer_type&& obs)
    : total_step_(tstep), step_count_(0), time_(0.), system_(std::move(sys)),
      ff_(std::move(ff)), integrator_(std::move(integr)), observer_(std::move(obs))
    {}
    ~MDSimulator() override = default;

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
    system_type     system_;
    forcefield_type ff_;
    integrator_type integrator_;
    observer_type   observer_;
};

template<typename traitsT>
inline void MDSimulator<traitsT>::initialize()
{
    this->integrator_.initialize(this->system_);
    this->ff_.initialize(this->system_, integrator_->delta_t());
    observer_.output(this->system_);
    return;
}

template<typename traitsT>
inline bool MDSimulator<traitsT>::step()
{
    this->time_ = integrator_.step(this->time_, system_, ff_);
    if(observer_.is_output_time())
    {
        observer_.output(this->time_, this->system_)
    }
    ++step_count_;
    return step_count <= total_step_;
}

template<typename traitsT>
inline void MDSimulator<traitsT>::finalize()
{
    observer_.output(this->system_);
    return;
}

} // mjolnir
#endif /* MJOLNIR_MD_SIMULATOR */
