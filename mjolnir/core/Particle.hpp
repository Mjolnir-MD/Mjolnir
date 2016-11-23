#ifndef MJOLNIR_PARTICLE
#define MJOLNIR_PARTICLE
#include <mjolnir/util/scalar_type_of.hpp>

namespace mjolnir
{

template<typename coordT>
struct Particle
{
    scalar_type<coordT> mass;
    coordT position;
    coordT velocity;
    coordT force;
};

} // mjolnir
#endif /* MJOLNIR_PARTICLE */
