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

template<typename realT, typename coordT>
constexpr inline Particle<realT, coordT>
make_particle(const realT  mass, const coordT& pos,
              const coordT& vel, const coordT& f)
{
    return Particle<realT, coordT>{mass, pos, vel, f};
}


} // mjolnir
#endif /* MJOLNIR_PARTICLE */
