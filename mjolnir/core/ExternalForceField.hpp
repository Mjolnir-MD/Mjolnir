#ifndef MJOLNIR_EXTERNAL_FORCE_FIELD
#define MJOLNIR_EXTERNAL_FORCE_FIELD
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/util/logger.hpp>
#include <utility>
#include <vector>
#include <array>
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class ExternalForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type         real_type;
    typedef typename traits_type::coordinate_type   coordinate_type;
    typedef System<traits_type>                     system_type;
    typedef ExternalForceInteractionBase<traitsT>   interaction_type;
    typedef std::unique_ptr<interaction_type>       interaction_ptr;
    typedef std::vector<interaction_ptr>            container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    ExternalForceField()  = default;
    ~ExternalForceField() = default;
    ExternalForceField(ExternalForceField const&) = delete;
    ExternalForceField(ExternalForceField&&)      = default;
    ExternalForceField& operator=(ExternalForceField const&) = delete;
    ExternalForceField& operator=(ExternalForceField&&)      = default;

    void emplace(interaction_ptr&& interaction)
    {
        interactions_.push_back(std::move(interaction));
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

} // mjolnir
#endif// MJOLNIR_EXTERNAL_FORCE_FIELD
