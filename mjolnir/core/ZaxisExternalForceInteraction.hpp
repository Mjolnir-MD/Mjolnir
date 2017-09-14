#ifndef MJOLNIR_GLOBAL_EXTERNAL_INTERACTION
#define MJOLNIR_GLOBAL_EXTERNAL_INTERACTION
#include "GlobalInteractionBase.hpp"
#include <memory>

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename partitionT>
class ZaxisExternalForceInteraction final : public GlobalInteractionBase<traitsT>
{
  public:
    
    using traits_type             = traitsT;
    using potential_type          = potentialT;
    using partition_type          = partitionT;
    using base_type               = GlobalInteractionBase<traits_type>;
    using real_type               = typename base_type::real_type;
    using coordinate_type         = typename base_type::coordinate_type;
    using system_type             = typename base_type::system_type;
    using particle_type           = typename base_type::particle_type;
    using boundary_type           = typename base_type::boundary_type;
    
  public:
    ZaxisExternalForceInteraction()  = default;
    ~ZaxisExternalForceInteraction() = default;
    
    ZaxisExternalForceInteraction(potential_type&& pot, partition_type&& part)
	: potential_(std::move(pot)), partition_(std::move(part))
    {}
    
    void initialize(const system_type& sys, const real_type dt) override
    {
	partition_.initialize(sys);
	return ;	
    }
    
    void      calc_force(system_type&)              override;
    real_type calc_energy(const system_type&) const override;
    
  private:
    
    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT, typename potT, typename spaceT>
void ZaxisExternalForceInteraction<traitsT, potT, spaceT>::calc_force(
    system_type& sys)
{
    partition_.update(sys);
    for(std::size_t i : partition_.obj_indices())
    {
	const real_type force = -potential_.derivative(i, sys[i].position.at(2));
	sys[i].force.at(2) += force;
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename ZaxisExternalForceInteraction<traitsT, potT, spaceT>::real_type 
ZaxisExternalForceInteraction<traitsT, potT, spaceT>::calc_energy(
    const system_type& sys) const
{
    real_type e = 0.0;
    for(std::size_t i = 0; i < sys.size(); ++i){
	e += potential_.potential(i, sys[i].position.at(2));
    }
    return e;
}

} //mjolnir
#endif /* MJOLNIR_GLOBAL_EXTERNAL_INTERACTION */
