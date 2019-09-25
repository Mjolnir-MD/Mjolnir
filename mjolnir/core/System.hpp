#ifndef MJOLNIR_CORE_SYSTEM_HPP
#define MJOLNIR_CORE_SYSTEM_HPP
#include <mjolnir/core/Unit.hpp>
#include <mjolnir/core/Particle.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <mjolnir/util/logger.hpp>
#include <vector>
#include <map>
#include <cassert>

namespace mjolnir
{

template<typename traitsT>
class System
{
  public:
    using traits_type     = traitsT;
    using string_type     = std::string;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using boundary_type   = typename traits_type::boundary_type;
    using topology_type   = Topology;
    using attribute_type  = std::map<std::string, real_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;

    using particle_type            = Particle<real_type, coordinate_type>;
    using particle_view_type       = ParticleView<real_type, coordinate_type>;
    using particle_const_view_type = ParticleConstView<real_type, coordinate_type>;

    using real_container_type          = std::vector<real_type>;
    using coordinate_container_type    = std::vector<coordinate_type>;

  public:

    System(const std::size_t num_particles, const boundary_type& bound)
        : velocity_initialized_(false), boundary_(bound),
          topology_(num_particles),  attributes_(),
          num_particles_(num_particles), masses_   (num_particles),
          rmasses_      (num_particles), positions_(num_particles),
          velocities_   (num_particles), forces_   (num_particles)
    {}
    ~System() = default;
    System(const System&) = default;
    System(System&&)      = default;
    System& operator=(const System&) = default;
    System& operator=(System&&)      = default;

    void initialize(rng_type& rng)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(this->velocity_initialized_)
        {
            MJOLNIR_LOG_NOTICE(
                "velocity is already given, nothing to initialize in System");
            return ;
        }
        if(!this->has_attribute("temperature"))
        {
            throw std::runtime_error("[error] to generate velocity, "
                    "system.attributes.temperature is required.");
        }

        const real_type kB    = physics::constants<real_type>::kB();
        const real_type T_ref = this->attribute("temperature");

        MJOLNIR_LOG_NOTICE("generating velocity with T = ", T_ref, "...");

        // generate Maxwell-Boltzmann distribution
        const real_type kBT = kB * T_ref;
        for(std::size_t i=0; i<this->size(); ++i)
        {
            const auto vel_coef = std::sqrt(kBT / this->mass(i));
            math::X(this->velocity(i)) = rng.gaussian(0, vel_coef);
            math::Y(this->velocity(i)) = rng.gaussian(0, vel_coef);
            math::Z(this->velocity(i)) = rng.gaussian(0, vel_coef);
        }
        MJOLNIR_LOG_NOTICE("done.");
        return;
    }

    coordinate_type adjust_direction(coordinate_type dr) const noexcept
    {return boundary_.adjust_direction(dr);}
    coordinate_type  adjust_position(coordinate_type dr) const noexcept
    {return boundary_.adjust_position(dr);}

    std::size_t size() const noexcept {return num_particles_;}

    particle_view_type operator[](std::size_t i) noexcept
    {
        return particle_view_type{
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_[i],
            this->topology_.name_of (i, std::nothrow),
            this->topology_.group_of(i, std::nothrow)
        };
    }
    particle_const_view_type operator[](std::size_t i) const noexcept
    {
        return particle_const_view_type{
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_[i],
            this->topology_.name_of (i, std::nothrow),
            this->topology_.group_of(i, std::nothrow)
        };
    }
    particle_view_type at(std::size_t i)
    {
        return particle_view_type{
            masses_.at(i),    rmasses_.at(i),
            positions_.at(i), velocities_.at(i), forces_.at(i),
            topology_.name_of(i), topology_.group_of(i)
        };
    }
    particle_const_view_type at(std::size_t i) const
    {
        return particle_const_view_type{
            masses_.at(i),    rmasses_.at(i),
            positions_.at(i), velocities_.at(i), forces_.at(i),
            topology_.name_of(i), topology_.group_of(i)
        };
    }

    real_type  mass (std::size_t i) const noexcept {return masses_[i];}
    real_type& mass (std::size_t i)       noexcept {return masses_[i];}
    real_type  rmass(std::size_t i) const noexcept {return rmasses_[i];}
    real_type& rmass(std::size_t i)       noexcept {return rmasses_[i];}

    coordinate_type const& position(std::size_t i) const noexcept {return positions_[i];}
    coordinate_type&       position(std::size_t i)       noexcept {return positions_[i];}
    coordinate_type const& velocity(std::size_t i) const noexcept {return velocities_[i];}
    coordinate_type&       velocity(std::size_t i)       noexcept {return velocities_[i];}
    coordinate_type const& force   (std::size_t i) const noexcept {return forces_[i];}
    coordinate_type&       force   (std::size_t i)       noexcept {return forces_[i];}

    string_type const& name (std::size_t i) const noexcept {return topology_.name_of(i, std::nothrow);}
    string_type&       name (std::size_t i)       noexcept {return topology_.name_of(i, std::nothrow);}
    string_type const& group(std::size_t i) const noexcept {return topology_.group_of(i, std::nothrow);}
    string_type&       group(std::size_t i)       noexcept {return topology_.group_of(i, std::nothrow);}

    boundary_type&       boundary()       noexcept {return boundary_;}
    boundary_type const& boundary() const noexcept {return boundary_;}
    topology_type&       topology()       noexcept {return topology_;}
    topology_type const& topology() const noexcept {return topology_;}

    // system attributes like `reference temperature`, `ionic strength`, ...
    // assuming it will not be called so often.
    real_type  attribute(const std::string& key) const {return attributes_.at(key);}
    real_type& attribute(const std::string& key)       {return attributes_[key];}
    bool   has_attribute(const std::string& key) const {return attributes_.count(key) == 1;}

    bool  velocity_initialized() const noexcept {return velocity_initialized_;}
    bool& velocity_initialized()       noexcept {return velocity_initialized_;}

  private:

    bool           velocity_initialized_;
    boundary_type  boundary_;
    topology_type  topology_;
    attribute_type attributes_;

    std::size_t                  num_particles_;
    real_container_type          masses_;
    real_container_type          rmasses_; // r for reciprocal
    coordinate_container_type    positions_;
    coordinate_container_type    velocities_;
    coordinate_container_type    forces_;
    // names and groups are in Topology class
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class System<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class System<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class System<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class System<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_SYSTEM_HPP
