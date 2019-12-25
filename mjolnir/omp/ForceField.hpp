#ifndef MJOLNIR_OMP_FORCE_FIELD_HPP
#define MJOLNIR_OMP_FORCE_FIELD_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/core/ForceField.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
class ForceField<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type              = OpenMPSimulatorTraits<realT, boundaryT>;
    using real_type                = typename traits_type::real_type;
    using coordinate_type          = typename traits_type::coordinate_type;
    using system_type              = System<traits_type>;
    using local_forcefield_type    = LocalForceField<traits_type>;
    using global_forcefield_type   = GlobalForceField<traits_type>;
    using external_forcefield_type = ExternalForceField<traits_type>;

  public:

    ForceField(local_forcefield_type&& local, global_forcefield_type&& global)
        : local_(std::move(local)), global_(std::move(global))
    {}
    ForceField(local_forcefield_type&&    local,
               global_forcefield_type&&   global,
               external_forcefield_type&& external)
        : local_(std::move(local)), global_(std::move(global)),
          external_(std::move(external))
    {}
    ForceField()  = default;
    ~ForceField() = default;
    ForceField(const ForceField&) = default;
    ForceField(ForceField&&)      = default;
    ForceField& operator=(const ForceField&) = default;
    ForceField& operator=(ForceField&&)      = default;

    // this modify system::topology by using local interaction info.
    void initialize(system_type& sys)
    {
        // first, fetch current topology
        local_.write_topology(sys.topology());
        sys.topology().construct_molecules();

        // based on the topology, make exclusion list
        local_   .initialize(sys);
        global_  .initialize(sys);
        external_.initialize(sys);
    }

    // update parameters like temperature, ionic concentration, etc...
    void update(const system_type& sys)
    {
        local_   .update(sys);
        global_  .update(sys);
        external_.update(sys);
    }

    // update margin of neighbor list
    void reduce_margin(const real_type dmargin, const system_type& sys)
    {
        local_   .reduce_margin(dmargin, sys);
        global_  .reduce_margin(dmargin, sys);
        external_.reduce_margin(dmargin, sys);
    }
    void scale_margin(const real_type scale, const system_type& sys)
    {
        local_   .scale_margin(scale, sys);
        global_  .scale_margin(scale, sys);
        external_.scale_margin(scale, sys);
        return;
    }

    void calc_force(system_type& sys) const noexcept
    {
#pragma omp parallel
        {
            // calc_force uses `nowait` to speedup. To do that, it needs to be
            // inside a parallel region. So here, we need to wrap `calc_force`
            // by a `parallel` region.
            local_   .calc_force(sys);
            global_  .calc_force(sys);
            external_.calc_force(sys);
        }
        // merge thread-local forces into master
        sys.merge_forces();
    }
    real_type calc_energy(const system_type& sys) const noexcept
    {
        return local_.calc_energy(sys) + global_.calc_energy(sys) +
            external_.calc_energy(sys);
    }

    std::vector<std::string> list_energy_name() const
    {
        auto retval = local_.list_energy();
        auto glo    = global_.list_energy();
        auto ext    = external_.list_energy();

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(std::make_move_iterator(glo.begin()),
                  std::make_move_iterator(glo.end()),
                  std::back_inserter(retval));

        std::copy(std::make_move_iterator(ext.begin()),
                  std::make_move_iterator(ext.end()),
                  std::back_inserter(retval));

        return retval;
    }
    std::vector<real_type> dump_energy(const system_type& sys) const
    {
        auto retval = local_.dump_energy(sys);
        auto glo    = global_.dump_energy(sys);
        auto ext    = external_.dump_energy(sys);

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
        std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
        return retval;
    }

    local_forcefield_type    const& local()    const noexcept {return local_;}
    local_forcefield_type    &      local()          noexcept {return local_;}
    global_forcefield_type   const& global()   const noexcept {return global_;}
    global_forcefield_type   &      global()         noexcept {return global_;}
    external_forcefield_type const& external() const noexcept {return external_;}
    external_forcefield_type &      external()       noexcept {return external_;}

  private:

    local_forcefield_type    local_;
    global_forcefield_type   global_;
    external_forcefield_type external_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ForceField<OpenMPSimulatorTraits<double, UnlimitedBoundary>>;
extern template class ForceField<OpenMPSimulatorTraits<float,  UnlimitedBoundary>>;
extern template class ForceField<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ForceField<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
