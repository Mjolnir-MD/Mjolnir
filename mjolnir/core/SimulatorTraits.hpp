#ifndef MJOLNIR_DEFAULT_TRAITS
#define MJOLNIR_DEFAULT_TRAITS
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/Particle.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
struct SimulatorTraitsBase
{
    typedef realT                                 real_type;
    typedef Vector<real_type, 3>                  coordinate_type;
    typedef coordinate_type                       position_type;
    typedef Particle<coordinate_type>             particle_type;
    typedef boundaryT<real_type, coordinate_type> boundary_type;
};

} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
