#ifndef MJOLNIR_CORE_SYSTEM_HPP
#define MJOLNIR_CORE_SYSTEM_HPP
#include <mjolnir/core/Unit.hpp>
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
    using matrix33_type   = typename traits_type::matrix33_type;
    using boundary_type   = typename traits_type::boundary_type;
    using topology_type   = Topology;
    using attribute_type  = std::map<std::string, real_type>;
    using rng_type        = RandomNumberGenerator<traits_type>;

    using real_container_type          = std::vector<real_type>;
    using coordinate_container_type    = std::vector<coordinate_type>;
    using string_container_type        = std::vector<std::string>;

  public:

    System(const std::size_t num_particles, const boundary_type& bound)
        : velocity_initialized_(false), force_initialized_(false),
          boundary_(bound), attributes_(), virial_(0,0,0, 0,0,0, 0,0,0),
          num_particles_(num_particles), masses_   (num_particles),
          rmasses_      (num_particles), positions_(num_particles),
          velocities_   (num_particles), forces_   (num_particles),
          names_        (num_particles), groups_   (num_particles)
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

        // make all the particles inside the boundary
        for(auto& p : this->positions_)
        {
            p = this->boundary_.adjust_position(p);
        }

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

    coordinate_type adjust_direction(coordinate_type from, coordinate_type to) const noexcept
    {
        return boundary_.adjust_direction(from, to);
    }
    coordinate_type  adjust_position(coordinate_type dr) const noexcept
    {
        return boundary_.adjust_position(dr);
    }
    coordinate_type transpose(coordinate_type tgt, const coordinate_type& ref) const noexcept
    {
        return boundary_.transpose(tgt, ref);
    }

    std::size_t size() const noexcept {return num_particles_;}

    // When parallelizing a code, forces are often calculated separately in
    // several computational units, like cores, nodes, gpu devices, etc. To
    // make it consistent, we may need to do something with forces calculated
    // separately. Those functions are provided for such an specialized
    // situation. Here, for the normal case, we do not need to do anything.
    //     Before calling `preprocess_forces()`, (a part of) forces may already
    // be calculated. So this function should NOT break the forces that is
    // already written in the `force(i)`.
    //     After calling `postprocess_forces()`, the result of `force(i)` always
    // represents the "force" of a particle at that time point. I mean, we can
    // consider the "force" is equivalent to the force that is calculated by
    // single core.
    void preprocess_forces()  noexcept {/* do nothing */}
    void postprocess_forces() noexcept {/* do nothing */}

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

    string_type const& name (std::size_t i) const noexcept {return names_[i];}
    string_type&       name (std::size_t i)       noexcept {return names_[i];}
    string_type const& group(std::size_t i) const noexcept {return groups_[i];}
    string_type&       group(std::size_t i)       noexcept {return groups_[i];}

    matrix33_type&       virial()       noexcept {return virial_;}
    matrix33_type const& virial() const noexcept {return virial_;}

    boundary_type&       boundary()       noexcept {return boundary_;}
    boundary_type const& boundary() const noexcept {return boundary_;}

    // system attributes like `reference temperature`, `ionic strength`, ...
    // assuming it will not be called so often.
    real_type  attribute(const std::string& key) const {return attributes_.at(key);}
    real_type& attribute(const std::string& key)       {return attributes_[key];}
    bool   has_attribute(const std::string& key) const {return attributes_.count(key) == 1;}
    attribute_type const& attributes() const noexcept {return attributes_;}

    bool  velocity_initialized() const noexcept {return velocity_initialized_;}
    bool& velocity_initialized()       noexcept {return velocity_initialized_;}
    bool  force_initialized()    const noexcept {return force_initialized_;}
    bool& force_initialized()          noexcept {return force_initialized_;}

    coordinate_container_type const& forces() const noexcept {return forces_;}
    coordinate_container_type&       forces()       noexcept {return forces_;}

  private:

    bool           velocity_initialized_, force_initialized_;
    boundary_type  boundary_;
    attribute_type attributes_;
    matrix33_type  virial_;

    std::size_t                  num_particles_;
    real_container_type          masses_;
    real_container_type          rmasses_; // r for reciprocal
    coordinate_container_type    positions_;
    coordinate_container_type    velocities_;
    coordinate_container_type    forces_;
    string_container_type        names_;
    string_container_type        groups_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own System<OpenMP> to avoid data race.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class System<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class System<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class System<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class System<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_SYSTEM_HPP
