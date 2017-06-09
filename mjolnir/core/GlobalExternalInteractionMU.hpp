#ifndef MJOLNIR_GLOBAL_EXTERNAL_INTERACTION
#define MJOLNIR_GLOBAL_EXTERNAL_INTERACTION
#include "GlobalInteractionBase.hpp"

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename partition>
class GlobalExternalInteraction : public GlobalInteractionBase<traitsT>
{
  public:
    
    using traits_type             = traitsT;
    using potential_type          = potentialT;
    using partition_type = partitionT;
    using base_type               = GlobalInteractionBase<traits_type>;
    using real_type               = base_type::real_type;
    using coordinate_type         = base_type::coordinate_type;
    using system_type             = base_type::system_type;
    using particle_type           = base_type::particle_type;
    using boundary_type           = base_type::boundary_type;
    
  public:
    GlobalExternalInteraction()  = default;
    ~GlobalExternalInteraction() = default;
    
    GlobalExternalInteraction(potential_type&& pot, partition_type&& part)
	: potential_(std::move(pot)),spatial_partition_(std::move(part))
    {}
    
    void initialize(const system_type&, const real_type) override
    {
	return ;
	
    }
    
    void calc_force(system_type&)                          override;
    void calc_energy(const particle_container_type&) const override;
    
  private:
    
    potential_type potential_;
    spatial_partition_type spatial_partition_;
};

template<typename traitsT, typename potT, typename spaceT>
void GlobalExternalInteraction<traitsT, potT, spaceT>::calc_force(
    system_type& sys)
{
    for(std::size_t i = 0; i < sys.size(); ++i){
	const real_type force = potential_.derivative(sys[i].position.at(2));
	sys[i].force.at(2) += f;
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename GlobalExternalInteraction<traitsT, potT, spaceT>::calc_energy
(system_type& sys) const
{
    real_type e = 0.0;
    for(std::size_t i = 0; i < sys.size(); ++i){
	e += potential_.potential(sys[i].position.at(2));
    }
    return e;
}

}  
#endif /* MJOLNIR_GLOBAL_EXTERNAL_INTERACTION */
