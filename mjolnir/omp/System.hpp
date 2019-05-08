#ifndef MJOLNIR_OMP_SYSTEM_HPP
#define MJOLNIR_OMP_SYSTEM_HPP
#include <mjolnir/util/aligned_allocator.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
class System<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using boundary_type   = typename traits_type::boundary_type;
    using topology_type   = Topology;
    using attribute_type  = std::map<std::string, real_type>;

    using particle_type            = Particle<real_type, coordinate_type>;
    using static_string_type       = typename particle_type::static_string_type;
    using particle_view_type       = ParticleView<real_type, coordinate_type>;
    using particle_const_view_type = ParticleConstView<real_type, coordinate_type>;

    using real_container_type          = std::vector<real_type>;
    using coordinate_container_type    = std::vector<coordinate_type>;
    using static_string_container_type = std::vector<static_string_type>;

    static constexpr std::size_t cache_alignment = 64;

    template<typename T>
    using cache_aligned_allocator = aligned_allocator<T, cache_alignment>;

  public:

    System(const std::size_t num_particles, const boundary_type& bound)
        : boundary_(bound), topology_(num_particles),
          attributes_(), num_particles_(num_particles),
          masses_   (num_particles), rmasses_   (num_particles),
          positions_(num_particles), velocities_(num_particles),
          forces_master_(num_particles),
          forces_threads_(omp_get_max_threads(),
              coordinate_container_type(num_particles),
              cache_aligned_allocator<coordinate_container_type>{}),
          names_    (num_particles), groups_    (num_particles)
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
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_master_[i],
            names_[i],     groups_[i]
        };
    }
    particle_const_view_type operator[](std::size_t i) const noexcept
    {
        return particle_const_view_type{
            masses_[i],    rmasses_[i],
            positions_[i], velocities_[i], forces_master_[i],
            names_[i],     groups_[i]
        };
    }
    particle_view_type at(std::size_t i)
    {
        return particle_view_type{
            masses_.at(i),    rmasses_.at(i),
            positions_.at(i), velocities_.at(i), forces_master_.at(i),
            names_.at(i),     groups_.at(i)
        };
    }
    particle_const_view_type at(std::size_t i) const
    {
        return particle_const_view_type{
            masses_.at(i),    rmasses_.at(i),
            positions_.at(i), velocities_.at(i), forces_master_.at(i),
            names_.at(i),     groups_.at(i)
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
    coordinate_type const& force   (std::size_t i) const noexcept {return forces_master_[i];}
    coordinate_type&       force   (std::size_t i)       noexcept {return forces_master_[i];}

    coordinate_type const&
    force_thread(std::size_t thread_num, std::size_t particle_id) const noexcept
    {
        return forces_threads_[thread_num][particle_id];
    }
    coordinate_type&
    force_thread(std::size_t thread_num, std::size_t particle_id)       noexcept
    {
        return forces_threads_[thread_num][particle_id];
    }

    static_string_type const& name(std::size_t i) const noexcept {return names_[i];}
    static_string_type&       name(std::size_t i)       noexcept {return names_[i];}
    static_string_type const& group(std::size_t i) const noexcept {return groups_[i];}
    static_string_type&       group(std::size_t i)       noexcept {return groups_[i];}

    boundary_type&       boundary()       noexcept {return boundary_;}
    boundary_type const& boundary() const noexcept {return boundary_;}
    topology_type&       topology()       noexcept {return topology_;}
    topology_type const& topology() const noexcept {return topology_;}

    // system attributes like `reference temperature`, `ionic strength`, ...
    // assuming it will not be called so often.
    real_type  attribute(const std::string& key) const {return attributes_.at(key);}
    real_type& attribute(const std::string& key)       {return attributes_[key];}
    bool   has_attribute(const std::string& key) const {return attributes_.count(key) == 1;}

  private:

    boundary_type  boundary_;
    topology_type  topology_;
    attribute_type attributes_;

    std::size_t                  num_particles_;
    real_container_type          masses_;
    real_container_type          rmasses_; // r for reciprocal
    coordinate_container_type    positions_;
    coordinate_container_type    velocities_;
    coordinate_container_type    forces_master_;
    std::vector<coordinate_container_type,
                cache_aligned_allocator<coordinate_container_type>
        > forces_threads_;
    static_string_container_type names_;
    static_string_container_type groups_;
};

} // mjolnir
#endif// MJOLNIR_SYSTEM_HPP
