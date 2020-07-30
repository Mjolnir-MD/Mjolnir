#ifndef MJOLNIR_OMP_BOND_LENGTH_GO_CONTACT_INTERACTION_HPP
#define MJOLNIR_OMP_BOND_LENGTH_GO_CONTACT_INTERACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/omp/System.hpp>
#include <mjolnir/forcefield/local/BondLengthGoContactInteraction.hpp>

namespace mjolnir
{

template<typename realT, template<typename, typename> class boundaryT>
class BondLengthInteraction<
    OpenMPSimulatorTraits<realT, boundaryT>,
    GoContactPotential<realT>
    > final : public LocalInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type          = OpenMPSimulatorTraits<realT, boundaryT>;
    using potential_type       = GoContactPotential<realT>;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using indices_type         = std::array<std::size_t, 2>;
    using potential_index_pair = std::pair<indices_type, potential_type>;
    using container_type       = std::vector<potential_index_pair>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

  public:

    BondLengthInteraction() = default;
    ~BondLengthInteraction() override {}

    BondLengthInteraction(const connection_kind_type kind,
                          const container_type& pot)
        : kind_(kind), potentials_(pot)
    {}
    BondLengthInteraction(const connection_kind_type kind,
                          container_type&& pot)
        : kind_(kind), potentials_(std::move(pot))
    {}

    void calc_force(system_type& sys) const noexcept override
    {
#pragma omp parallel for
        for(std::size_t i=0; i<this->potentials_.size(); ++i)
        {
            const auto& idxp = this->potentials_[i];

            const std::size_t idx0 = idxp.first[0];
            const std::size_t idx1 = idxp.first[1];
            const auto&       pot  = idxp.second;

            const auto dpos =
                sys.adjust_direction(sys.position(idx1) - sys.position(idx0));

            const real_type len2  = math::length_sq(dpos);
            if(pot.cutoff() * pot.cutoff() <= len2)
            {
                continue;
            }

            const real_type r2     = real_type(1) / len2;
            const real_type v0r_2  = pot.v0() * pot.v0() * r2;
            const real_type v0r_6  = v0r_2 * v0r_2 * v0r_2;
            const real_type v0r_10 = v0r_6 * v0r_2 * v0r_2;
            const real_type v0r_12 = v0r_10 * v0r_2;

            const auto coef = -60 * pot.k() * r2 * (v0r_10 - v0r_12);
            const auto f    = coef * dpos;

            const std::size_t thread_id = omp_get_thread_num();
            sys.force_thread(thread_id, idx0) -= f;
            sys.force_thread(thread_id, idx1) += f;
        }
        return;
    }
    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E = 0;
#pragma omp parallel for reduction(+:E)
        for(std::size_t i=0; i<this->potentials_.size(); ++i)
        {
            const auto& idxp = this->potentials_[i];
            E += idxp.second.potential(math::length(sys.adjust_direction(
                    sys.position(idxp.first[1]) - sys.position(idxp.first[0]))));
        }
        return E;
    }
    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        real_type E = 0;
#pragma omp parallel for reduction(+:E)
        for(std::size_t i=0; i<this->potentials_.size(); ++i)
        {
            const auto& idxp = this->potentials_[i];

            const std::size_t idx0 = idxp.first[0];
            const std::size_t idx1 = idxp.first[1];
            const auto&       pot  = idxp.second;

            const auto dpos =
                sys.adjust_direction(sys.position(idx1) - sys.position(idx0));

            const real_type len2  = math::length_sq(dpos);
            if(pot.cutoff() * pot.cutoff() <= len2)
            {
                continue;
            }

            const real_type r2     = real_type(1) / len2;
            const real_type v0r_2  = pot.v0() * pot.v0() * r2;
            const real_type v0r_6  = v0r_2 * v0r_2 * v0r_2;
            const real_type v0r_10 = v0r_6 * v0r_2 * v0r_2;
            const real_type v0r_12 = v0r_10 * v0r_2;

            E += pot.k() * (5 * v0r_12 - 6 * v0r_10);

            const auto coef = -60 * pot.k() * r2 * (v0r_10 - v0r_12);
            const auto f    = coef * dpos;

            const std::size_t thread_id = omp_get_thread_num();
            sys.force_thread(thread_id, idx0) -= f;
            sys.force_thread(thread_id, idx1) += f;
        }
        return E;
    }


    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of bonds = ", potentials_.size());
        for(auto& potential : potentials_)
        {
            potential.second.initialize(sys);
        }
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials_)
        {
            item.second.update(sys);
        }
    }

    // do nothing. this is used to reduce margin of neighbor list, and added
    // to this class for the consistency.
    void reduce_margin(const real_type, const system_type&) override {return;}
    void  scale_margin(const real_type, const system_type&) override {return;}

    std::string name() const override {return "BondLengthGoContact"_s;}

    void write_topology(topology_type& topol) const override
    {
        if(this->kind_.empty() || this->kind_ == "none") {return;}

        for(const auto& idxp : this->potentials_)
        {
            const auto i = idxp.first[0];
            const auto j = idxp.first[1];
            topol.add_connection(i, j, this->kind_);
        }
        return;
    }

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new BondLengthInteraction(kind_, container_type(potentials_));
    }

  private:
    connection_kind_type kind_;
    container_type potentials_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{

// go-contact
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, UnlimitedBoundary       >, GoContactPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  UnlimitedBoundary       >, GoContactPotential<float> >;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<double, CuboidalPeriodicBoundary>, GoContactPotential<double>>;
extern template class BondLengthInteraction<OpenMPSimulatorTraits<float,  CuboidalPeriodicBoundary>, GoContactPotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_INTERACTION_BOND_LENGTH_GO_CONTACT_INTERACTION_HPP */
