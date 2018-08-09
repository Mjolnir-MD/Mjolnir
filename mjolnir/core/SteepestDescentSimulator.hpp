#ifndef MJOLNIR_STEEPEST_DESCENT_SIMULATOR
#define MJOLNIR_STEEPEST_DESCENT_SIMULATOR
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/core/Observer.hpp>
#include <limits>

namespace mjolnir
{

template<typename traitsT>
class SteepestDescentSimulator final : public SimulatorBase
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type>     system_type;
    typedef ForceField<traits_type> forcefield_type;
    typedef Observer<traits_type>   observer_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    SteepestDescentSimulator(const real_type h, const real_type threshold,
            const std::size_t step_limit, const std::size_t save_step,
            system_type&& sys, forcefield_type&& ff, observer_type&& obs)
    : h_(h), threshold_(threshold),
      step_limit_(step_limit), step_count_(0), save_step_(save_step),
      system_(std::move(sys)), ff_(std::move(ff)), observer_(std::move(obs))
    {}
    ~SteepestDescentSimulator() override = default;

    void initialize() override;
    bool step()       override;
    void finalize()   override;

    real_type calc_energy() const {return this->ff_.calc_energy(this->system_);}

    system_type&       system()       noexcept {return system_;}
    system_type const& system() const noexcept {return system_;}

    ForceField<traitsT>&       forcefields()       noexcept {return ff_;}
    ForceField<traitsT> const& forcefields() const noexcept {return ff_;}

  protected:
    real_type       h_;
    real_type       threshold_;
    std::size_t     step_limit_;
    std::size_t     step_count_;
    std::size_t     save_step_;
    system_type     system_;
    forcefield_type ff_;
    observer_type   observer_;
};

template<typename traitsT>
inline void SteepestDescentSimulator<traitsT>::initialize()
{
    this->ff_.initialize(this->system_);

    this->observer_.initialize(this->system_, this->ff_);
    this->observer_.output(0, this->system_, this->ff_);
    return;
}

template<typename traitsT>
inline bool SteepestDescentSimulator<traitsT>::step()
{
    // calculate negative derivatives (-dV/dr)
    this->ff_.calc_force(this->system_);

    real_type max_disp2 = 0.0; // to update cell list and check the convergence
    real_type max_diff  = 0.0; // to check the convergence
    for(std::size_t i=0; i<this->system_.size(); ++i)
    {
        const coordinate_type disp = this->h_ * this->system_[i].force;

        max_diff = std::max(max_diff, std::abs(this->system_[i].force[0]));
        max_diff = std::max(max_diff, std::abs(this->system_[i].force[1]));
        max_diff = std::max(max_diff, std::abs(this->system_[i].force[2]));

        max_disp2 = std::max(max_disp2, length_sq(disp));
        system_[i].position = system_.adjust_position(system_[i].position + disp);
        system_[i].force    = coordinate_type(0, 0, 0);
    }

    if(max_diff < this->threshold_)
    {
        this->observer_.output(static_cast<real_type>(this->step_count_),
                               this->system_, this->ff_);
        return false; // converged. stop the simulation!
    }
    this->system_.largest_displacement() = std::sqrt(max_disp2);

    if(step_count_ % save_step_ == 0)
    {
        this->observer_.output(static_cast<real_type>(this->step_count_),
                               this->system_, this->ff_);
    }
    ++step_count_;
    return this->step_count_ < this->step_limit_;
}

template<typename traitsT>
inline void SteepestDescentSimulator<traitsT>::finalize()
{
    this->observer_.output(static_cast<real_type>(this->step_count_),
                           this->system_, this->ff_);
    return;
}

} // mjolnir
#endif /* MJOLNIR_STEEPEST_DESCENT_SIMULATOR */
