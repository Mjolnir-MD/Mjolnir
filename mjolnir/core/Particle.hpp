#ifndef MJOLNIR_PARTICLE
#define MJOLNIR_PARTICLE
#include <mjolnir/util/static_string.hpp>

namespace mjolnir
{

template<typename realT, typename coordT>
struct Particle
{
    typedef realT            real_type;
    typedef coordT           coordinate_type;
    typedef static_string<8> static_string_type;
    real_type          mass;
    real_type          rmass; // r is for reciprocal
    coordinate_type    position;
    coordinate_type    velocity;
    coordinate_type    force;
    static_string_type name;
};

// To use SoA data structure
// class System {
//     std::vector<Real>  masses;
//     std::vector<Real>  rmasses;
//     std::vector<Coord> positions;
//     std::vector<Coord> velocities;
//     std::vector<Coord> forces;
// }
// ParticleView and ConstView provides AoS-like accessibility.

template<typename realT, typename coordT>
class ParticleView
{
  public:
    typedef realT  real_type;
    typedef coordT coordinate_type;
    typedef Particle<real_type, coordinate_type>       particle_type;
    typedef typename particle_type::static_string_type static_string_type;

    ParticleView()  = delete; // it makes all the references undefined.
    ~ParticleView() = default;
    ParticleView(const ParticleView&) = default;
    ParticleView(ParticleView&&)      = default;
    ParticleView& operator=(const ParticleView&) = default;
    ParticleView& operator=(ParticleView&&)      = default;

    ParticleView(real_type& m, real_type& rm,
        coordinate_type& p, coordinate_type& v,
        coordinate_type& f, static_string_type& n) noexcept
        : mass(m), rmass(rm), position(p), velocity(v), force(f), name(n)
    {}

    ParticleView& operator=(const Particle<real_type, coordinate_type>& p)
    {
        mass     = p.mass;
        rmass    = p.rmass;
        position = p.position;
        velocity = p.velocity;
        force    = p.force;
        name     = p.name;
        return *this;
    }

    real_type&          mass;
    real_type&          rmass;
    coordinate_type&    position;
    coordinate_type&    velocity;
    coordinate_type&    force;
    static_string_type& name;
};

template<typename realT, typename coordT>
class ParticleConstView
{
  public:
    typedef realT  real_type;
    typedef coordT coordinate_type;
    typedef Particle<real_type, coordinate_type>       particle_type;
    typedef typename particle_type::static_string_type static_string_type;

    ParticleConstView()  = delete; // it makes all the references undefined.
    ~ParticleConstView() = default;
    ParticleConstView(const ParticleConstView&) = default;
    ParticleConstView(ParticleConstView&&)      = default;
    ParticleConstView& operator=(const ParticleConstView&) = default;
    ParticleConstView& operator=(ParticleConstView&&)      = default;

    ParticleConstView(const real_type& m, const real_type& rm,
        const coordinate_type& p, const coordinate_type& v,
        const coordinate_type& f, const static_string_type& n) noexcept
        : mass(m), rmass(rm), position(p), velocity(v), force(f), name(n)
    {}

    real_type          const& mass;
    real_type          const& rmass;
    coordinate_type    const& position;
    coordinate_type    const& velocity;
    coordinate_type    const& force;
    static_string_type const& name;
};

} // mjolnir
#endif /* MJOLNIR_PARTICLE */
