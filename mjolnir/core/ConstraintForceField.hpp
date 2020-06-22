#ifndef MJOLNIR_CORE_CONSTRAINT_FORCE_FIELD_HPP
#define MJOLNIR_CORE_CONSTRAINT_FORCE_FIELD_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/util/logger.hpp>

namespace mjolnir
{

template<typename traitsT>
class ConstraintForceField
{

// Constraint class for g-BAOAB Langevin integrator developed by the following papers
// Leimkuhler B., Matthews C., Proc. R. Soc. A Math. Phys. Eng. Sci. (2016)
  public:
    using traits_type          = traitsT; 
    using system_type          = System<traits_type>;
    using real_type            = typename traits_type::real_type;
    using topology_type        = typename system_type::topology_type;
    using connection_kind_type = typename topology_type::connection_kind_type;

    using indices_v0_pair_type = std::pair<std::array<std::size_t, 2>, real_type>;
    using container_type       = std::vector<indices_v0_pair_type>;

  public:

    ConstraintForceField()  = default;
    ConstraintForceField(const connection_kind_type kind,
                         const container_type&& constraints,
                         const std::size_t      max_iteration,
                         const real_type        tolerance)
        : max_iteration_(max_iteration), tolerance_(tolerance), kind_(kind), 
          constraints_(std::move(constraints))
    {}
    ~ConstraintForceField() = default;
    ConstraintForceField(ConstraintForceField&&)            = default;
    ConstraintForceField& operator=(ConstraintForceField&&) = default;
    ConstraintForceField(const ConstraintForceField& other) = default;

    void write_topology(typename system_type::topology_type& topol)
    {
        if(this->kind_.empty() || this->kind_ == "none") {return;}

        for(auto& idxc : this->constraints_)
        {
            const auto i = idxc.first[0];
            const auto j = idxc.first[1];
            topol.add_connection(i, j, this->kind_);
        }
        return;
    }

    std::size_t           max_iteration() const noexcept {return max_iteration_;}
    real_type             tolerance()     const noexcept {return tolerance_;}
    const container_type& constraints()   const noexcept {return constraints_;}

    bool                  empty()         const noexcept {return constraints_.empty();}

  private:
    std::size_t          max_iteration_;
    real_type            tolerance_;
    connection_kind_type kind_;

    container_type         constraints_;
};

#ifdef MJOLNIR_SEPARATE_BUILD 
extern template class ConstraintForceField<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class ConstraintForceField<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class ConstraintForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class ConstraintForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif /* MJOLNIR_GLOBAL_FORCE_FIELD */
