#ifndef MJOLNIR_SIMULATOR
#define MJOLNIR_SIMULATOR
#include "Integrator.hpp"

namespace mjolnir
{

template<typename traitsT>
class Simulator
{
  public:
    typedef typename traitsT::real_type real_type;
    typedef typename traitsT::time_type time_type;

    Simulator(ParticleContainer<traitsT>&& pcon,
              ForceField<traitsT>&& ff,
              std::unique_ptr<Integrator<traitsT>>&& integr)
        : time_(0), pcon_(std::forward<ParticleContainer<traitsT>>(pcon)),
          ff_(std::forward<ForceField<traitsT>>(ff)),
          integrator_(std::forward<std::unique_ptr<Integrator<traitsT>>>(integr))
    {}
    ~Simulator() = default;

    void initialize();
    void step();
    void step_until(const time_type t);
    real_type calc_energy() const;

    ParticleContainer<traitsT>&       particles()       {return pcon_;}
    ParticleContainer<traitsT> const& particles() const {return pcon_;}

    ForceField<traitsT>&       forcefields()       {return ff_;}
    ForceField<traitsT> const& forcefields() const {return ff_;}

    time_type& time()       {return time_;}
    time_type  time() const {return time_;}

  protected:
    time_type time_;
    ParticleContainer<traitsT>           pcon_;
    ForceField<traitsT>                  ff_;
    std::unique_ptr<Integrator<traitsT>> integrator_;
};

template<typename traitsT>
inline void Simulator<traitsT>::initialize()
{
    this->integrator_->initialize(this->pcon_);
    this->ff_.initialize(this->pcon_, integrator_->delta_t());
    return;
}

template<typename traitsT>
inline void Simulator<traitsT>::step()
{
    this->time_ = integrator_->step(this->time_, pcon_, ff_);
    return;
}

template<typename traitsT>
inline void Simulator<traitsT>::step_until(const time_type t)
{
    while(this->time_ < t)
    {
        this->time_ = integrator_->step(this->time_, pcon_, ff_);
    }
    return ;
}

template<typename traitsT>
inline typename Simulator<traitsT>::real_type
Simulator<traitsT>::calc_energy() const
{
    return this->ff_.calc_energy(this->pcon_);
}


} // mjolnir
#endif /* MJOLNIR_SIMULATOR */
