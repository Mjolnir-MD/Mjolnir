#ifndef MJOLNIR_PARTICLE
#define MJOLNIR_PARTICLE

namespace mjolnir
{

template<typename realT, typename coordT>
struct Particle
{
    typedef realT  real_type;
    typedef coordT coordinate_type;
    real_type mass;
    coordinate_type position;
    coordinate_type velocity;
    coordinate_type force;
};

} // mjolnir
#endif /* MJOLNIR_PARTICLE */
