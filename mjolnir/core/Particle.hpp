#ifndef MJOLNIR_PARTICLE
#define MJOLNIR_PARTICLE

namespace mjolnir
{

template<typename realT, typename coordT>
struct Particle
{
    typedef realT  real_type;
    typedef coordT coordinate_type;
    real_type       mass;
    coordinate_type position;
    coordinate_type velocity;
    coordinate_type force;
};

// To use SoA data structure
// class System {
//     std::vector<Real>  masses;
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

    ParticleView()  = delete; // it makes all the references undefined.
    ~ParticleView() = default;
    ParticleView(const ParticleView&) = default;
    ParticleView(ParticleView&&)      = default;
    ParticleView& operator=(const ParticleView&) = default;
    ParticleView& operator=(ParticleView&&)      = default;

    ParticleView(real_type&       m, coordinate_type& p,
                 coordinate_type& v, coordinate_type& f) noexcept
        : mass(m), position(p), velocity(v), force(f)
    {}

    ParticleView& operator=(const Particle<real_type, coordinate_type>& p)
    {
        mass     = p.mass;
        position = p.position;
        velocity = p.velocity;
        force    = p.force;
        return *this;
    }

    real_type&       mass;
    coordinate_type& position;
    coordinate_type& velocity;
    coordinate_type& force;
};

template<typename realT, typename coordT>
class ParticleConstView
{
  public:
    typedef realT  real_type;
    typedef coordT coordinate_type;

    ParticleConstView()  = delete; // it makes all the references undefined.
    ~ParticleConstView() = default;
    ParticleConstView(const ParticleConstView&) = default;
    ParticleConstView(ParticleConstView&&)      = default;
    ParticleConstView& operator=(const ParticleConstView&) = default;
    ParticleConstView& operator=(ParticleConstView&&)      = default;

    ParticleConstView(const real_type&       m, const coordinate_type& p,
                      const coordinate_type& v, const coordinate_type& f) noexcept
        : mass(m), position(p), velocity(v), force(f)
    {}

    real_type       const& mass;
    coordinate_type const& position;
    coordinate_type const& velocity;
    coordinate_type const& force;
};

} // mjolnir
#endif /* MJOLNIR_PARTICLE */
