#ifndef MJOLNIR_SYSTEM
#define MJOLNIR_SYSTEM
#include "StructureTopology.hpp"
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
    typedef typename traits_type::particle_type     particle_type;
    typedef StructureTopology                       topology_type;
    typedef std::map<std::string, real_type>        attribute_type;
    typedef std::vector<particle_type>              container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

    System()  = default;
    ~System() = default;

    System(container_type&& pcon, boundary_type&& bound)
        : boundary_(bound), particles_(pcon)
    {
        this->topology_.resize(this->particles_.size());
    }

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

    real_type& max_speed()       noexcept {return max_speed_;}
    real_type  max_speed() const noexcept {return max_speed_;}

  private:
    real_type      max_speed_;
    boundary_type  boundary_;
    container_type particles_;
    topology_type  topology_;
    attribute_type attributes_;
};

} // mjolnir
#endif// MJOLNIR_SYSTEM
