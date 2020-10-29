#ifndef MJOLNIR_OMP_BOND_LENGTH_INTERACTION_HPP
#define MJOLNIR_OMP_BOND_LENGTH_INTERACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/forcefield/local/BondLengthInteraction.hpp>

namespace mjolnir
{

/*! @brief calculate energy and force of Bond length type local interaction */
template<typename realT, template<typename, typename> class boundaryT,
         typename potentialT>
class BondLengthInteraction<OpenMPSimulatorTraits<realT, boundaryT>, potentialT>
    final : public LocalInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type          = OpenMPSimulatorTraits<realT, boundaryT>;
    using potential_type       = potentialT;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using indices_type         = std::array<std::size_t, 2>;
    using potential_index_pair = std::pair<indices_type, potentialT>;
    using container_type       = std::vector<potential_index_pair>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

  public:

    BondLengthInteraction(const connection_kind_type kind,
                          const container_type& pot)
        : kind_(kind), potentials(pot)
    {}
    BondLengthInteraction(const connection_kind_type kind,
                          container_type&& pot)
        : kind_(kind), potentials(std::move(pot))
    {}
    ~BondLengthInteraction() override {}

    void      calc_force (system_type& sys)       const noexcept override
    {
#pragma omp parallel for
        for(std::size_t i=0; i<this->potentials.size(); ++i)
        {
            const auto& idxp = this->potentials[i];

            const std::size_t idx0 = idxp.first[0];
            const std::size_t idx1 = idxp.first[1];

            const auto dpos =
                sys.adjust_direction(sys.position(idx0), sys.position(idx1));

            const real_type len2 = math::length_sq(dpos); // l^2
            const real_type rlen = math::rsqrt(len2);     // 1/l
            const real_type force = -1 * idxp.second.derivative(len2 * rlen);
            // here, L^2 * (1 / L) = L.

            const coordinate_type f = dpos * (force * rlen);
            sys.force_thread(omp_get_thread_num(), idx0) -= f;
            sys.force_thread(omp_get_thread_num(), idx1) += f;
        }
        return;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E = 0.;
#pragma omp parallel for reduction(+:E)
        for(std::size_t i=0; i<this->potentials.size(); ++i)
        {
            const auto& idxp = this->potentials[i];
            E += idxp.second.potential(math::length(sys.adjust_direction(
                    sys.position(idxp.first[0]), sys.position(idxp.first[1]))));
        }
        return E;
    }

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        real_type E = 0.;
#pragma omp parallel for reduction(+:E)
        for(std::size_t i=0; i<this->potentials.size(); ++i)
        {
            const auto& idxp = this->potentials[i];

            const std::size_t idx0 = idxp.first[0];
            const std::size_t idx1 = idxp.first[1];

            const auto dpos =
                sys.adjust_direction(sys.position(idx0), sys.position(idx1));

            const real_type len2 = math::length_sq(dpos); // l^2
            const real_type rlen = math::rsqrt(len2);     // 1/l
            const real_type  len = len2 * rlen;
            const real_type force = -1 * idxp.second.derivative(len);
            E += idxp.second.potential(len);

            const coordinate_type f = dpos * (force * rlen);
            sys.force_thread(omp_get_thread_num(), idx0) -= f;
            sys.force_thread(omp_get_thread_num(), idx1) += f;
        }
        return E;
    }

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("With OpenMP: potential = ", potential_type::name(),
                         ", number of bonds = ", potentials.size());
        for(auto& potential : potentials)
        {
            potential.second.initialize(sys);
        }
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials)
        {
            item.second.update(sys);
        }
    }

    // do nothing. this is used to reduce margin of neighbor list, and added
    // to this class for the consistency.
    void reduce_margin(const real_type, const system_type&) override {return;}
    void  scale_margin(const real_type, const system_type&) override {return;}

    std::string name() const override
    {return "BondLength:"_s + potential_type::name();}

    void write_topology(topology_type& topol) const override
    {
        if(this->kind_.empty() || this->kind_ == "none") {return;}

        for(const auto& idxp : this->potentials)
        {
            const auto i = idxp.first[0];
            const auto j = idxp.first[1];
            topol.add_connection(i, j, this->kind_);
        }
        return;
    }

    base_type* clone() const override
    {
        return new BondLengthInteraction(kind_, container_type(potentials));
    }

  private:
    connection_kind_type kind_;
    container_type potentials;
};

} // mjolnir

#include <mjolnir/omp/BondLengthGoContactInteraction.hpp>

#ifdef MJOLNIR_SEPARATE_BUILD
// explicitly specialize BondLengthInteraction with LocalPotentials
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/forcefield/local/GoContactPotential.hpp>
#include <mjolnir/forcefield/local/GaussianPotential.hpp>

namespace mjolnir
{

// harmonic
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, HarmonicPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, HarmonicPotential<float> >;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// Specialization for go-contact is defined in BondLengthGoContactInteraction.hpp.

// gaussian
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
