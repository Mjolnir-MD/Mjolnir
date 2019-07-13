#ifndef MJOLNIR_CORE_OBSERVER_CONTAINER_HPP
#define MJOLNIR_CORE_OBSERVER_CONTAINER_HPP
#include <mjolnir/util/io.hpp>
#include <mjolnir/util/progress_bar.hpp>
#include <mjolnir/core/ObserverBase.hpp>
#include <vector>

// This class manages several different XXXObservers to output
// positions, energies, topologies, and others.
//
// Also, this class outputs progress bar into stdout if the flag is on.

namespace mjolnir
{

template<typename traitsT>
class ObserverContainer
{
  public:
    using observer_base_type = ObserverBase<traitsT>;
    using observer_base_ptr  = std::unique_ptr<observer_base_type>;
    using real_type          = typename observer_base_type::real_type;
    using coordinate_type    = typename observer_base_type::coordinate_type;
    using system_type        = typename observer_base_type::system_type;
    using forcefield_type    = typename observer_base_type::forcefield_type;
    using progress_bar_type  = progress_bar</*width of bar = */50>;

  public:

    explicit ObserverContainer(bool output_progress = false)
        : output_progress_(output_progress && io::detail::isatty(std::cerr))
    {}
    ~ObserverContainer() = default;

    ObserverContainer(const ObserverContainer&) = default;
    ObserverContainer(ObserverContainer&&)      = default;
    ObserverContainer& operator=(const ObserverContainer&) = default;
    ObserverContainer& operator=(ObserverContainer&&)      = default;

    void initialize(const std::size_t total_step, const real_type dt,
                    const system_type& sys, const forcefield_type& ff)
    {
        for(const auto& obs : observers_)
        {
            obs->initialize(total_step, dt, sys, ff);
        }

        this->progress_bar_.reset(total_step); // set total_step as 100%.
    }
    void output(const std::size_t step, const real_type dt,
                const system_type& sys, const forcefield_type& ff)
    {
        for(const auto& obs : this->observers_)
        {
            obs->output(step, dt, sys, ff);
        }

        // this branching might be wiped out by introducing another parameter
        // to SimulatorTraits, but I don't think the cost of the implementation
        // is larger than the benefit on the runtime efficiency.
        if(this->output_progress_)
        {
            std::cerr << this->progress_bar_.format(step);
        }
    }
    void finalize(const std::size_t total_step, const real_type dt,
                  const system_type& sys, const forcefield_type& ff)
    {
        for(const auto& obs : this->observers_)
        {
            obs->finalize(total_step, dt, sys, ff);
        }

        if(this->output_progress_)
        {
            // In order to re-write progress bar in each step, it does not print
            // end-line after the progress bar. But if finalize is called, we
            // can `finalize` the progress bar.
            std::cerr << std::endl;
        }
    }

    // assign one another XXXObserver
    void push_back(observer_base_ptr&& obs)
    {
        this->observers_.push_back(std::move(obs));
    }

    // mainly for testing purpose.
    std::vector<observer_base_ptr> const& observers() const noexcept
    {
        return this->observers_;
    }

  private:
    std::vector<observer_base_ptr> observers_;
    progress_bar_type              progress_bar_;
    bool                           output_progress_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ObserverContainer<SimulatorTraits<double, UnlimitedBoundary>       >;
extern template class ObserverContainer<SimulatorTraits<float,  UnlimitedBoundary>       >;
extern template class ObserverContainer<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ObserverContainer<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_CORE_OBSERVER_HPP
