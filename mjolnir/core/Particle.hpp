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

template<typename realT, typename coordT>
constexpr inline Particle<realT, coordT>
make_particle(const realT  mass, const coordT& pos,
              const coordT& vel, const coordT& f)
{
    return Particle<realT, coordT>{mass, pos, vel, f};
}


} // mjolnir
#endif /* MJOLNIR_PARTICLE */
