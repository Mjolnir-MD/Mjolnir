#ifndef MJOLNIR_PARTICLE
#define MJOLNIR_PARTICLE

namespace mjolnir
{

template<typename realT, typename coordT>
struct Particle
{
    typedef realT  real_type;
    typedef coordT coordinate_type;
    realT  mass;
    coordT position;
    coordT velocity;
    coordT force;
};

} // mjolnir
#endif /* MJOLNIR_PARTICLE */
