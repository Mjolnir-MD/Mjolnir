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
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef Topology                              topology_type;
    typedef std::map<std::string, real_type>      attribute_type;

    typedef Particle<real_type, coordinate_type>          particle_type;
    typedef ParticleView<real_type, coordinate_type>      particle_view_type;
    typedef ParticleConstView<real_type, coordinate_type> particle_const_view_type;

  public:

    System(const std::size_t num_particles, const boundary_type& bound)
        : largest_disp_(0.0), boundary_(bound), topology_(num_particles),
          attributes_(), num_particles_(num_particles),
          masses_   (num_particles), rmasses_   (num_particles),
          positions_(num_particles), velocities_(num_particles),
          forces_   (num_particles)
    {}
    ~System() = default;

    coordinate_type adjust_direction(coordinate_type dr) const noexcept
    {return boundary_.adjust_direction(dr);}
    coordinate_type  adjust_position(coordinate_type dr) const noexcept
    {return boundary_.adjust_position(dr);}

    std::size_t size() const noexcept {return num_particles_;}

    particle_view_type operator[](std::size_t i) noexcept
    {
        return particle_view_type{
            masses_[i], rmasses_[i], positions_[i], velocities_[i], forces_[i]
        };
    }
    particle_const_view_type operator[](std::size_t i) const noexcept
    {
        return particle_const_view_type{
            masses_[i], rmasses_[i], positions_[i], velocities_[i], forces_[i]
        };
    }
    particle_view_type at(std::size_t i)
    {
        return particle_view_type{masses_.at(i), rmasses_.at(i),
            positions_.at(i), velocities_.at(i), forces_.at(i)
        };
    }
    particle_const_view_type at(std::size_t i) const
    {
        return particle_const_view_type{masses_.at(i), rmasses_.at(i),
            positions_.at(i), velocities_.at(i), forces_.at(i)
        };
    }

    std::vector<real_type>&             masses()           noexcept {return masses_;}
    std::vector<real_type> const&       masses()     const noexcept {return masses_;}
    std::vector<real_type>&             rmasses()          noexcept {return masses_;}
    std::vector<real_type> const&       rmasses()    const noexcept {return masses_;}
    std::vector<coordinate_type>&       positions()        noexcept {return positions_;}
    std::vector<coordinate_type> const& positions()  const noexcept {return positions_;}
    std::vector<coordinate_type>&       velocities()       noexcept {return velocities_;}
    std::vector<coordinate_type> const& velocities() const noexcept {return velocities_;}
    std::vector<coordinate_type>&       forces()           noexcept {return forces_;}
    std::vector<coordinate_type> const& forces()     const noexcept {return forces_;}

    // to implement these, we need a ZipIterator. but is it really needed?
//     iterator       begin()        noexcept {return particles_.begin();}
//     iterator       end()          noexcept {return particles_.end();}
//     const_iterator begin()  const noexcept {return particles_.cbegin();}
//     const_iterator end()    const noexcept {return particles_.cend();}
//     const_iterator cbegin() const noexcept {return particles_.cbegin();}
//     const_iterator cend()   const noexcept {return particles_.cend();}

    boundary_type&       boundary()       noexcept {return boundary_;}
    boundary_type const& boundary() const noexcept {return boundary_;}
    topology_type&       topology()       noexcept {return topology_;}
    topology_type const& topology() const noexcept {return topology_;}

    // system attributes like `reference temperature`, `ionic strength`, ...
    // assuming it will not be called so often.
    real_type  attribute(const std::string& key) const {return attributes_.at(key);}
    real_type& attribute(const std::string& key)       {return attributes_[key];}
    bool   has_attribute(const std::string& key) const {return attributes_.count(key) == 1;}

    // store largest displacement in a step to decrease margin of verlet list.
    real_type& largest_displacement()       noexcept {return largest_disp_;}
    real_type  largest_displacement() const noexcept {return largest_disp_;}

  private:
    real_type      largest_disp_;
    boundary_type  boundary_;
    topology_type  topology_;
    attribute_type attributes_;

    std::size_t num_particles_;
    std::vector<real_type>       masses_;
    std::vector<real_type>       rmasses_; // r for reciprocal
    std::vector<coordinate_type> positions_;
    std::vector<coordinate_type> velocities_;
    std::vector<coordinate_type> forces_;
};

} // mjolnir
#endif// MJOLNIR_SYSTEM
