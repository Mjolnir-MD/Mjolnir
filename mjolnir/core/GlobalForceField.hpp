#ifndef MJOLNIR_CORE_GLOBAL_FORCE_FIELD_HPP
#define MJOLNIR_CORE_GLOBAL_FORCE_FIELD_HPP
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <array>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class GlobalForceField
{
  public:
    using traits_type      = traitsT;
    using real_type        = typename traits_type::real_type;
    using coordinate_type  = typename traits_type::coordinate_type;
    using boundary_type    = typename traits_type::boundary_type;
    using system_type      = System<traits_type>;
    using interaction_base = GlobalInteractionBase<traitsT>;
    using interaction_ptr  = std::unique_ptr<interaction_base>;
    using container_type   = std::vector<interaction_ptr>;
    using iterator         = typename container_type::iterator;
    using const_iterator   = typename container_type::const_iterator;

  public:
    GlobalForceField() = default;
    ~GlobalForceField() = default;
    GlobalForceField(const GlobalForceField&) = delete;
    GlobalForceField(GlobalForceField&&)      = default;
    GlobalForceField& operator=(const GlobalForceField&) = delete;
    GlobalForceField& operator=(GlobalForceField&&)      = default;

    void emplace(interaction_ptr&& inter)
    {
        interactions_.emplace_back(std::move(inter));
    }

    void initialize(const system_type& sys)
    {
        for(auto& item : this->interactions_)
        {
            item->initialize(sys);
        }
    }

    // to re-calculate parameters like temperature, ionic concentration, etc...
    void update(const system_type& sys)
    {
        for(auto& item : this->interactions_)
        {
            item->update(sys);
        }
    }

    // to reduce margin of neighbor list, and re-construct the list if needed
    void update_margin(const real_type dmargin, const system_type& sys)
    {
        for(auto& item : this->interactions_)
        {
            item->update_margin(dmargin, sys);
        }
    }

    void calc_force(system_type& sys) const noexcept
    {
        for(const auto& item : this->interactions_)
        {
            item->calc_force(sys);
        }
        return;
    }
    real_type calc_energy(const system_type& sys) const noexcept
    {
        real_type energy = 0.;
        for(const auto& item : this->interactions_)
        {
            energy += item->calc_energy(sys);
        }
        return energy;
    }

    // basically, it is called only once at the begenning of a simulation.
    // this function do a lot of stuff, such as memory allocation, but it does
    // not affect runtime efficiency so much.
    std::vector<std::string> list_energy() const
    {
        std::vector<std::string> retval;
        for(const auto& i : interactions_)
        {
            retval.push_back(i->name());
        }
        return retval;
    }

    std::vector<real_type> dump_energy(const system_type& sys) const
    {
        std::vector<real_type> retval;
        for(const auto& i : interactions_)
        {
            retval.push_back(i->calc_energy(sys));
        }
        return retval;
    }

    bool           empty()  const noexcept {return interactions_.empty();}
    std::size_t    size()   const noexcept {return interactions_.size();}
    iterator       begin()        noexcept {return interactions_.begin();}
    iterator       end()          noexcept {return interactions_.end();}
    const_iterator begin()  const noexcept {return interactions_.begin();}
    const_iterator end()    const noexcept {return interactions_.end();}
    const_iterator cbegin() const noexcept {return interactions_.begin();}
    const_iterator cend()   const noexcept {return interactions_.end();}

  private:

    container_type interactions_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class GlobalForceField<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class GlobalForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class GlobalForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class GlobalForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_GLOBAL_FORCE_FIELD */
