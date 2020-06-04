#ifndef MJOLNIR_OMP_SYSTEM_HPP
#define MJOLNIR_OMP_SYSTEM_HPP
#include <mjolnir/util/aligned_allocator.hpp>
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/RandomNumberGenerator.hpp>
#include <mjolnir/core/System.hpp>

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
    using rng_type        = RandomNumberGenerator<traits_type>;

    using string_type              = std::string;
    using particle_type            = Particle<real_type, coordinate_type>;
    using particle_view_type       = ParticleView<real_type, coordinate_type>;
    using particle_const_view_type = ParticleConstView<real_type, coordinate_type>;

    using real_container_type          = std::vector<real_type>;
    using coordinate_container_type    = std::vector<coordinate_type>;
    using string_container_type        = std::vector<std::string>;

    static constexpr std::size_t cache_alignment = 64;

    template<typename T>
    using cache_aligned_allocator = aligned_allocator<T, cache_alignment>;

  public:

    System(const std::size_t num_particles, const boundary_type& bound)
        : velocity_initialized_(false), boundary_(bound),
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
        assert(this->has_attribute("temperature"));

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

    void merge_forces() noexcept
    {
#pragma omp parallel for
        for(std::size_t i=0; i<this->size(); ++i)
        {
            // sumup all the thread-local forces into the main force.
            // XXX: keep the values in the main force! Some pre-omp interactions
            //      may write its force in it. Clearing main force is a task
            //      for Integrators.
            for(std::size_t thread_id=0, max_threads=omp_get_max_threads();
                    thread_id < max_threads; ++thread_id)
            {
                this->force(i) += this->force_thread(thread_id, i);
                this->force_thread(thread_id, i) =
                    math::make_coordinate<coordinate_type>(0, 0, 0);
            }
        }
        return;
    }

    string_type const& name (std::size_t i) const noexcept {return names_[i];}
    string_type&       name (std::size_t i)       noexcept {return names_[i];}
    string_type const& group(std::size_t i) const noexcept {return groups_[i];}
    string_type&       group(std::size_t i)       noexcept {return groups_[i];}

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

  private:

    bool           velocity_initialized_;
    boundary_type  boundary_;
    attribute_type attributes_;

    std::size_t                  num_particles_;
    real_container_type          masses_;
    real_container_type          rmasses_; // r for reciprocal
    coordinate_container_type    positions_;
    coordinate_container_type    velocities_;
    coordinate_container_type    forces_master_;
    // thread-local forces
    std::vector<coordinate_container_type,
                cache_aligned_allocator<coordinate_container_type>
        > forces_threads_;
    string_container_type        names_;
    string_container_type        groups_;

    // names and groups are in Topology class
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class System<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
extern template class System<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
extern template class System<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class System<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_SYSTEM_HPP
