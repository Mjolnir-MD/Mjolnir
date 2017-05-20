#ifndef MJOLNIR_MD_SIMULATOR
#define MJOLNIR_MD_SIMULATOR
#include "SimulatorBase.hpp"
#include "System.hpp"

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
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    MDSimulator(System<traitsT>&& sys, ForceField<traitsT>&& ff,
                integrator_type&& integr)
        : time_(0), system_(std::move(sys)), ff_(std::move(ff)),
          integrator_(std::move(integr))
    {}
    ~MDSimulator() override = default;

    void initialize() override;
    void step()       override;
    void finalize()   override {return;}

    real_type calc_energy() const;

    system_type&       particles()       {return system_;}
    system_type const& particles() const {return system_;}

    ForceField<traitsT>&       forcefields()       {return ff_;}
    ForceField<traitsT> const& forcefields() const {return ff_;}

    real_type& time()       {return time_;}
    real_type  time() const {return time_;}

  protected:
    real_type       time_;
    system_type     system_;
    forcefield_type ff_;
    integrator_type integrator_;
};

template<typename traitsT>
inline void MDSimulator<traitsT>::initialize()
{
    this->integrator_->initialize(this->pcon_);
    this->ff_.initialize(this->pcon_, integrator_->delta_t());
    return;
}

template<typename traitsT>
inline void MDSimulator<traitsT>::step()
{
    this->time_ = integrator_->step(this->time_, pcon_, ff_);
    return;
}

template<typename traitsT>
inline typename MDSimulator<traitsT>::real_type
MDSimulator<traitsT>::calc_energy() const
{
    return this->ff_.calc_energy(this->pcon_);
}


} // mjolnir
#endif /* MJOLNIR_MD_SIMULATOR */
