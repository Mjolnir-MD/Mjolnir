#ifndef MJOLNIR_CORE_PARTICLE_HPP
#define MJOLNIR_CORE_PARTICLE_HPP
#include <string>

namespace mjolnir
{

template<typename realT, typename coordT>
struct Particle
{
    using real_type       = realT;
    using coordinate_type = coordT;
    using string_type     = std::string;

    real_type       mass;
    real_type       rmass; // r is for reciprocal
    coordinate_type position;
    coordinate_type velocity;
    coordinate_type force;
    string_type     name;
    string_type     group;
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
    using particle_type   = Particle<realT, coordT>;
    using real_type       = typename particle_type::real_type;
    using coordinate_type = typename particle_type::coordinate_type;
    using string_type     = std::string;

    ParticleView()  = delete; // it makes all the references undefined.
    ~ParticleView() = default;
    ParticleView(const ParticleView&) = default;
    ParticleView(ParticleView&&)      = default;
    ParticleView& operator=(const ParticleView&) = default;
    ParticleView& operator=(ParticleView&&)      = default;

    ParticleView(real_type& m, real_type&         rm,
        coordinate_type&    p, coordinate_type&    v, coordinate_type& f,
        string_type&        n, string_type&        g) noexcept
    : mass(m), rmass(rm), position(p), velocity(v), force(f), name(n), group(g)
    {}

    ParticleView& operator=(const Particle<real_type, coordinate_type>& p)
    {
        mass     = p.mass;
        rmass    = p.rmass;
        position = p.position;
        velocity = p.velocity;
        force    = p.force;
        name     = p.name;
        group    = p.group;
        return *this;
    }

    real_type&       mass;
    real_type&       rmass;
    coordinate_type& position;
    coordinate_type& velocity;
    coordinate_type& force;
    string_type&     name;
    string_type&     group;
};

template<typename realT, typename coordT>
class ParticleConstView
{
  public:
    using particle_type   = Particle<realT, coordT>;
    using real_type       = typename particle_type::real_type;
    using coordinate_type = typename particle_type::coordinate_type;
    using string_type     = std::string;

    ParticleConstView()  = delete; // it makes all the references undefined.
    ~ParticleConstView() = default;
    ParticleConstView(const ParticleConstView&) = default;
    ParticleConstView(ParticleConstView&&)      = default;
    ParticleConstView& operator=(const ParticleConstView&) = default;
    ParticleConstView& operator=(ParticleConstView&&)      = default;

    ParticleConstView(const real_type& m, const real_type& rm,
        const coordinate_type& p, const coordinate_type& v,
        const coordinate_type& f,
        const string_type& n, const string_type& g) noexcept
    : mass(m), rmass(rm), position(p), velocity(v), force(f), name(n), group(g)
    {}

    real_type       const& mass;
    real_type       const& rmass;
    coordinate_type const& position;
    coordinate_type const& velocity;
    coordinate_type const& force;
    string_type     const& name;
    string_type     const& group;
};

} // mjolnir
#endif /* MJOLNIR_PARTICLE */
