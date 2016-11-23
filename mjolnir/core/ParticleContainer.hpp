#ifndef MJOLNIR_PARTICLE_CONTAINER
#define MJOLNIR_PARTICLE_CONTAINER
#include "Particle.hpp"
#include <vector>

namespace mjolnir
{

template<typename traitsT>
class ParticleContainer
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::position_type position_type;
    typedef Particle<position_type> particle_type;
    typedef std::vector<particle_type> container_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:
    ParticleContainer(std::size_t number_of_particles)
        : particles_(number_of_particles)
    {}
    ~ParticleContainer() = default;

    particle_type &      operator[](std::size_t i)       {return particles_[i];}
    particle_type const& operator[](std::size_t i) const {return particles_[i];}
    particle_type &      at(std::size_t i)       {return particles_.at(i);}
    particle_type const& at(std::size_t i) const {return particles_.at(i);}

    iterator begin(){return particles_.begin();}
    iterator end()  {return particles_.end();}
    const_iterator cbegin() const {return particles_.cbegin();}
    const_iterator cend()   const {return particles_.cend();}

  private:
    container_type particles_;
};




} // mjolnir
#endif /* MJOLNIR_PARTICLE_CONTAINER */
