#ifndef MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#define MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#include "GlobalInteractionBase.hpp"
#include <memory>

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename partitionT>
class GlobalDistanceInteraction final : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef partitionT partition_type;
    typedef GlobalInteractionBase<traitsT> base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::particle_type   particle_type;
    typedef typename base_type::boundary_type   boundary_type;

  public:
    GlobalDistanceInteraction()  = default;
    ~GlobalDistanceInteraction() = default;

    GlobalDistanceInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}

    void initialize(const system_type& sys, const real_type dt) override
    {
        this->partition_.update(sys, dt);
    }

    void      calc_force (system_type&)             override;
    real_type calc_energy(const system_type&) const override;

  private:

    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT, typename potT, typename spaceT>
void GlobalDistanceInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys)
{
    potential_.update(sys);
    partition_.update(sys);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        for(auto j : this->partition_.partners(i))
        {
            const coordinate_type rij =
                sys.adjust_direction(sys[j].position - sys[i].position);
            const real_type       l = length(rij);
            const coordinate_type f = rij * (potential_.derivative(i, j, l) / l);
            sys[i].force += f;
            sys[j].force -= f;
        }
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename GlobalDistanceInteraction<traitsT, potT, spaceT>::real_type
GlobalDistanceInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        for(auto j : this->partition_.partners(i))
        {
            const real_type l = length(
                sys.adjust_direction(sys[j].position - sys[i].position));
            e += potential_.potential(i, j, l);
        }
    }
    return e;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
