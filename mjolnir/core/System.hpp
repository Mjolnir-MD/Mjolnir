#ifndef MJOLNIR_SYSTEM
#define MJOLNIR_SYSTEM
#include <mjolnir/core/Particle.hpp>
#include <mjolnir/core/Topology.hpp>
#include <vector>
#include <map>

namespace mjolnir
{

template<typename traitsT>
class System
{
  public:
    typedef traitsT  traits_type;
    typedef typename traits_type::real_type         real_type;
    typedef typename traits_type::coordinate_type   coordinate_type;
    typedef typename traits_type::boundary_type     boundary_type;
    typedef Particle<real_type, coordinate_type>    particle_type;
    typedef Topology                                topology_type;
    typedef std::map<std::string, real_type>        attribute_type;
    typedef std::vector<particle_type>              container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    System(const std::size_t num_particles, const boundary_type& bound)
        : boundary_(bound), particles_(num_particles)
    {
        this->topology_.resize(num_particles);
    }
    System(std::vector<particle_type>&& ps, const boundary_type& bound)
        : boundary_(bound), particles_(ps)
    {
        this->topology_.resize(particles_.size());
    }

    ~System() = default;

    coordinate_type adjust_direction(coordinate_type dr) const noexcept
    {return boundary_.adjust_direction(dr);}
    coordinate_type  adjust_position(coordinate_type dr) const noexcept
    {return boundary_.adjust_position(dr);}

    std::size_t size() const noexcept {return particles_.size();}

    particle_type &      operator[](std::size_t i)       noexcept {return particles_[i];}
    particle_type const& operator[](std::size_t i) const noexcept {return particles_[i];}
    particle_type &      at(std::size_t i)       {return particles_.at(i);}
    particle_type const& at(std::size_t i) const {return particles_.at(i);}

    iterator       begin()        noexcept {return particles_.begin();}
    iterator       end()          noexcept {return particles_.end();}
    const_iterator begin()  const noexcept {return particles_.cbegin();}
    const_iterator end()    const noexcept {return particles_.cend();}
    const_iterator cbegin() const noexcept {return particles_.cbegin();}
    const_iterator cend()   const noexcept {return particles_.cend();}

    boundary_type&       boundary()       noexcept {return boundary_;}
    boundary_type const& boundary() const noexcept {return boundary_;}
    topology_type&       topology()       noexcept {return topology_;}
    topology_type const& topology() const noexcept {return topology_;}

    // system attributes like `reference temperature`, `ionic strength`, ...
    // assuming it will not be called so often.
    real_type  attribute(const std::string& key) const {return attributes_.at(key);}
    real_type& attribute(const std::string& key)       {return attributes_[key];}

    // store largest displacement in a step to decrease margin of verlet list.
    real_type& largest_displacement()       noexcept {return largest_disp_;}
    real_type  largest_displacement() const noexcept {return largest_disp_;}

  private:
    real_type      largest_disp_;
    boundary_type  boundary_;
    container_type particles_;
    topology_type  topology_;
    attribute_type attributes_;
};

} // mjolnir
#endif// MJOLNIR_SYSTEM
