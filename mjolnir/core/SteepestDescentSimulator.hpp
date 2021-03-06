#ifndef MJOLNIR_CORE_STEEPEST_DESCENT_SIMULATOR_HPP
#define MJOLNIR_CORE_STEEPEST_DESCENT_SIMULATOR_HPP
#include <mjolnir/core/SimulatorBase.hpp>
#include <mjolnir/core/ObserverContainer.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/MsgPackSaver.hpp>
#include <mjolnir/math/math.hpp>
#include <limits>

namespace mjolnir
{

template<typename traitsT>
class SteepestDescentSimulator final : public SimulatorBase
{
  public:
    using traits_type        = traitsT;
    using real_type          = typename traits_type::real_type;
    using coordinate_type    = typename traits_type::coordinate_type;
    using system_type        = System<traits_type>;
    using forcefield_type    = std::unique_ptr<ForceFieldBase<traits_type>>;
    using observer_type      = ObserverContainer<traits_type>;
    using saver_type      = MsgPackSaver<traits_type>;

    SteepestDescentSimulator(const real_type h, const real_type threshold,
        const std::size_t step_limit, const std::size_t save_step,
        const std::size_t checkpoint_step,
        system_type&& sys, forcefield_type&& ff, observer_type&& obs)
    : h_(h), threshold_(threshold), step_limit_(step_limit), step_count_(0),
      save_step_(save_step), checkpoint_(checkpoint_step),
      system_(std::move(sys)), ff_(std::move(ff)), observers_(std::move(obs)),
      saver_(observers_.prefix())
    {}
    ~SteepestDescentSimulator() override {}

    void initialize() override;
    bool step()       override;
    void run()        override;
    void finalize()   override;

    system_type&       system()       noexcept {return system_;}
    system_type const& system() const noexcept {return system_;}

    forcefield_type&       forcefields()       noexcept {return ff_;}
    forcefield_type const& forcefields() const noexcept {return ff_;}

  protected:
    real_type       h_;
    real_type       threshold_;
    std::size_t     step_limit_;
    std::size_t     step_count_;
    std::size_t     save_step_;
    std::size_t     checkpoint_;
    system_type     system_;
    forcefield_type ff_;
    observer_type   observers_;
    saver_type      saver_;
};

template<typename traitsT>
inline void SteepestDescentSimulator<traitsT>::initialize()
{
    // XXX: Because this simulator does not use velocity,
    //      it does not initialize System.
    this->ff_->initialize(this->system_);

    // here, steepest_descent method has no physical `time`.
    // There is nothing we can except filling it with zero or something
    // that works as a marker.
    this->observers_.initialize(this->step_limit_, this->save_step_,
       /* there is no dt, so */ h_, this->system_, this->ff_);
    return;
}

template<typename traitsT>
inline bool SteepestDescentSimulator<traitsT>::step()
{
    if(step_count_ % save_step_ == 0)
    {
        this->observers_.output(this->step_count_, /* dt */ real_type(0.0),
                                this->system_, this->ff_);
    }
    if(step_count_ % checkpoint_ == 0)
    {
        saver_.save(this->system_);
    }

    // calculate negative derivatives (-dV/dr)
    this->ff_->calc_force(this->system_);

    real_type max_disp2 = 0.0; // to update cell list and check the convergence
    real_type max_diff  = 0.0; // to check the convergence
    for(std::size_t i=0; i<this->system_.size(); ++i)
    {
        const coordinate_type disp = this->h_ * this->system_.force(i);

        max_diff = std::max(max_diff, std::abs(math::X(this->system_.force(i))));
        max_diff = std::max(max_diff, std::abs(math::Y(this->system_.force(i))));
        max_diff = std::max(max_diff, std::abs(math::Z(this->system_.force(i))));

        max_disp2 = std::max(max_disp2, math::length_sq(disp));
        system_.position(i) = system_.adjust_position(system_.position(i) + disp);
        system_.force   (i) = math::make_coordinate<coordinate_type>(0, 0, 0);
    }

    if(max_diff < this->threshold_)
    {
        return false; // converged. stop the simulation!
    }

    // update neighbor list; reduce margin, reconstruct the list if needed
    this->ff_->reduce_margin(2 * std::sqrt(max_disp2), this->system_);

    ++step_count_;
    return this->step_count_ < this->step_limit_;
}

template<typename traitsT>
inline void SteepestDescentSimulator<traitsT>::run()
{
    while(this->step()){/* do nothing */;}
    return;
}

template<typename traitsT>
inline void SteepestDescentSimulator<traitsT>::finalize()
{
    this->observers_.output  (this->step_count_, /* dt */ real_type(0.0),
                              this->system_, this->ff_);
    this->observers_.finalize(this->step_limit_, /* dt */ real_type(0.0),
                              this->system_, this->ff_);
    this->saver_.save(this->system_);
    return;
}

#ifdef MJOLNIR_SEPARATE_BUILD
// BAOAB
extern template class SteepestDescentSimulator<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class SteepestDescentSimulator<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class SteepestDescentSimulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class SteepestDescentSimulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif // SEPARATE_BUILD

} // mjolnir
#endif /* MJOLNIR_STEEPEST_DESCENT_SIMULATOR */
