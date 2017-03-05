#ifndef MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#define MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#include "GlobalInteractionBase.hpp"
#include "SpatialPartition.hpp"
#include "BoundaryCondition.hpp"
#include <memory>

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename partitionT,
         typename boundaryT = UnlimitedBoundary<traitsT>>
class GlobalDistanceInteraction : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef partitionT spatial_partition_type;
    typedef boundaryT  boundary_type;
    typedef GlobalInteractionBase<traitsT> base_type;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::particle_container_type particle_container_type;

  public:
    GlobalDistanceInteraction()  = default;
    ~GlobalDistanceInteraction() = default;

    GlobalDistanceInteraction(potential_type&& pot,
                              spatial_partition_type&& space)
        : potential_(std::forward<potential_type>(pot)),
          spatial_partition_(std::forward<spatial_partition_type>(space))
    {}

    void
    initialize(const particle_container_type& pcon, const time_type dt) override;

    void
    calc_force(particle_container_type& pcon) override;

    real_type
    calc_energy(const particle_container_type& pcon) const override;

    void
    reset_parameter(const std::string& name, const real_type val) override;

  private:
    potential_type         potential_;
    spatial_partition_type spatial_partition_;
};

template<typename traitsT, typename potT, typename spaceT, typename boundaryT>
void GlobalDistanceInteraction<traitsT, potT, spaceT, boundaryT>::calc_force(
        particle_container_type& pcon)
{
    spatial_partition_.update(pcon);
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spatial_partition_type::index_list const& partners =
            spatial_partition_.partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = boundary_type::adjust_direction(
                    pcon[j].position - pcon[i].position);
            const real_type       l = length(rij);
            const coordinate_type f = rij * (potential_.derivative(i, j, l) / l);
            pcon[i].force += f;
            pcon[j].force -= f;
        }
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT, typename boundaryT>
typename GlobalDistanceInteraction<traitsT, potT, spaceT, boundaryT>::real_type
GlobalDistanceInteraction<traitsT, potT, spaceT, boundaryT>::calc_energy(
        const particle_container_type& pcon) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spatial_partition_type::index_list const& partners =
            spatial_partition_.partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = boundary_type::adjust_direction(
                    pcon[j].position - pcon[i].position);
            const real_type l = length(rij);
            e += potential_.potential(i, j, l);
        }
    }
    return e;
}

template<typename traitsT, typename potT, typename spaceT, typename boundaryT>
void GlobalDistanceInteraction<traitsT, potT, spaceT, boundaryT>::initialize(
        const particle_container_type& pcon, const time_type dt)
{
    this->spatial_partition_.update(pcon, dt);
    return ;
}

template<typename traitsT, typename potT, typename spaceT, typename boundaryT>
void GlobalDistanceInteraction<traitsT, potT, spaceT, boundaryT>::reset_parameter(
        const std::string& name, const real_type val)
{
    this->potential_.reset_parameter(name, val);
    return ;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
