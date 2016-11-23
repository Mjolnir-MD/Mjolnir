#ifndef MJOLNIR_DEFAULT_TRAITS
#define MJOLNIR_DEFAULT_TRAITS
#include <mjolnir/math/Vector.hpp>

namespace mjolnir
{

struct DefaultTraits
{
    typedef double               real_type;
    typedef real_type            time_type;
    typedef real_type            energy_type;
    typedef Vector<real_type, 3> vector_type;
    typedef vector_type          coordinate_type;
    typedef coordinate_type      position_type;
};

} // mjolnir
#endif /* MJOLNIR_DEFAULT_TRAITS */
