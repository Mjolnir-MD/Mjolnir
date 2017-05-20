#ifndef MJOLNIR_DEFAULT_TRAITS
#define MJOLNIR_DEFAULT_TRAITS
#include <mjolnir/math/Vector.hpp>

namespace mjolnir
{

template<typename realT, typename boundaryT>
struct SimulatorTraitsBase
{
    typedef realT                real_type;
    typedef Vector<real_type, 3> coordinate_type;
    typedef coordinate_type      position_type;
    typedef boundaryT            boundary_type;
};

} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
