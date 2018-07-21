#ifndef MJOLNIR_LOCAL_FORCE_FIELD
#define MJOLNIR_LOCAL_FORCE_FIELD
#include "LocalInteractionBase.hpp"
#include <mjolnir/util/logger.hpp>
#include <utility>
#include <vector>
#include <array>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class LocalForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type         real_type;
    typedef typename traits_type::coordinate_type   coordinate_type;
    typedef System<traits_type>                     system_type;
    typedef LocalInteractionBase<traitsT>           interaction_type;
    typedef std::unique_ptr<interaction_type>       interaction_ptr;
    typedef std::vector<interaction_ptr>            container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    LocalForceField()  = default;
    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

    void emplace(interaction_ptr&& interaction)
    {
        interactions_.emplace_back(std::move(interaction));
    }

    void initialize(const system_type& sys, const real_type dt)
    {
        for(auto& item : this->interactions_)
        {
            item->initialize(sys, dt);
        }
    }

    // to re-calculate parameters like temperature, ionic concentration, etc...
    void update(const system_type& sys)
    {
        for(auto& item : this->interactions_)
        {
            item->update(sys);
        }
        return;
    }

    // Topology is defined based on LocalForceField.
    void write_topology(typename system_type::topology_type& topol)
    {
        for(auto& item : this->interactions_)
        {
            item->write_topology(topol);
        }
        return;
    }

    void calc_force(system_type& sys) const
    {
        for(const auto& item : this->interactions_)
        {
            item->calc_force(sys);
        }
        return;
    }
    real_type calc_energy(const system_type& sys) const
    {
        real_type energy = 0.0;
        for(const auto& item : this->interactions_)
        {
            energy += item->calc_energy(sys);
        }
        return energy;
    }

    // TODO simplify
    std::string list_energy() const
    {
        std::string retval;
        for(const auto& i : interactions_)
        {
            retval += i->name();
            retval += ' ';
        }
        return retval;
    }

    std::string dump_energy(const system_type& sys) const
    {
        std::string retval;
        for(const auto& i : interactions_)
        {
            retval += std::to_string(i->calc_energy(sys));
            retval += ' ';
        }
        return retval;
    }

    iterator       begin()        noexcept {return interactions_.begin();}
    iterator       end()          noexcept {return interactions_.end();}
    const_iterator begin()  const noexcept {return interactions_.begin();}
    const_iterator end()    const noexcept {return interactions_.end();}
    const_iterator cbegin() const noexcept {return interactions_.begin();}
    const_iterator cend()   const noexcept {return interactions_.end();}

  private:

    container_type interactions_;
};

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
