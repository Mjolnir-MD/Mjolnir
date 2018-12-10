#ifndef MJOLNIR_GLOBAL_FORCE_FIELD
#define MJOLNIR_GLOBAL_FORCE_FIELD
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
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef GlobalInteractionBase<traitsT>    interaction_base;
    typedef std::unique_ptr<interaction_base> interaction_ptr;

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

    void calc_force(system_type& sys)
    {
        for(const auto& item : this->interactions_)
            item->calc_force(sys);
        return;
    }
    real_type calc_energy(const system_type& sys) const
    {
        real_type energy = 0.;
        for(const auto& item : this->interactions_)
            energy += item->calc_energy(sys);
        return energy;
    }

    // TODO simplify
    std::string list_energy() const
    {
        std::string retval;
        for(const auto& i : interactions_)
        {
            retval += ' ';
            retval += i->name();
        }
        return retval;
    }

    std::string dump_energy(const system_type& sys) const
    {
        std::ostringstream oss;
        for(const auto& i : interactions_)
        {
            oss << ' ' << std::setw(i->name().size()) << std::fixed
                << std::right << i->calc_energy(sys);
        }
        return oss.str();
    }

  private:

    std::vector<interaction_ptr> interactions_;
};

} // mjolnir
#endif /* MJOLNIR_GLOBAL_FORCE_FIELD */
