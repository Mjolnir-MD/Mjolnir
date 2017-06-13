#ifndef MJOLNIR_GLOBAL_EXTERNAL_INTERACTION
#define MJOLNIR_GLOBAL_EXTERNAL_INTERACTION
#include "GlobalInteractionBase.hpp"

namespace mjolnir
{

  template<typename traitsT, typename potentialT, typename partitionT,
	   typename boundaryT = UnlimitedBoundary<traitT>>
    class GlobalExternalInteraction : public GlobalInteractionBase<traitsT>
    {
    public:
    
      using traits_type = traitsT;
      using potential_type = potentialT;
      using spartial_partition_type = partitionT;
      using boundary_type = boudaryT;
      using base_type = GlobalInteractionBase<traits_type>;
      using time_type = base_type::time_type;
      using real_type = base_type::real_type;
      using coordinate_type = base_type::coordinate_type;
      using particle_container_type = base_type::particle_container_type;
    
    public:
      GlobalExternalInteraction()  = default;
      ~GlobalExternalInteraction() = default;

      GlobalExternalInteraction(potential_type&& pot,
				       spatial_partition_type&& space)
	: potential_(std::move<potential_type>(pot)),
	  spatial_partition_(std::move<spatial_partition_type>(space))
      {}

      void
      initialize(const particle_container_type& pcon, const time_type dt) override;
      void
      calc_force(particle_container_type& pcon) override;
      void
      calc_energy(const particle_container_type& pcon) override;
       
      
    private:
      potential_type potential_;
      spatial_partition_type spatial_partition_;
    };
  
  //TODO この先に関数の実装を書く
  template<typename traitT, typename potT, typename spaceT, typename boundaryT>
  void GlobalExternalInteraction<traitsT, potT, spaceT, boundaryT>::calc_force
  (particle_container_type& pcon)
  {
    spartial_partition_.update(pcon);//その時の系の状態に応じてpconを作りなおしている（要確認）
    typename spatial_partition_type::index_list const& partners =
      spatial_partition_.partners()//Prob ここで必要のない引数を指定してしまうことになりそう
    
  }

#endif /* MJOLNIR_GLOBAL_EXTERNAL_INTERACTION */
