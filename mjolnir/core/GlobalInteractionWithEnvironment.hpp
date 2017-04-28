#ifdef MJOLNIR_GLOBAL_INTERACTION_WITH_ENVIRONMENT
#define MJOLNIR_GLOBAL_INTERACTION_WITH_ENVIRONMENT
#include "GlobalInteractionBase.hpp"

namespace mjolnir
{

  template<typename traitT, typename potentialT, typename partitionT,
	   typename boundaryT = UnlimitedBoundary<traitT>>
    class GlobalInteractionWithEnvironment : public GlobalInteractionBase<traitT>
    {
    public:
      
    }
}
