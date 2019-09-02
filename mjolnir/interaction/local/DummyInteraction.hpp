#ifndef MJOLNIR_INTERACTION_DUMMY_INTERACTION_HPP
#define MJOLNIR_INTERACTION_DUMMY_INTERACTION_HPP
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/util/empty.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/string.hpp>

namespace mjolnir
{

// XXX DummyInteraction does not calculate force or energy. It is a dummy.
//
// Some potential has a strange topology depending on a property of a particle.
// This class is introduced to handle those special stuff in an easy way with
// the smallest overhead.
template<typename traitsT>
class DummyInteraction final : public LocalInteractionBase<traitsT>
{
  public:
    using traits_type          = traitsT;
    using potential_type       = empty_t; // XXX
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using indices_type         = std::array<std::size_t, 2>;
    using container_type       = std::vector<indices_type>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

  public:

    DummyInteraction(const connection_kind_type& kind, const container_type& cont)
        : kind_(kind), potentials_(cont)
    {}
    ~DummyInteraction() override {}

    void      calc_force (system_type&)       const noexcept override {return;}
    real_type calc_energy(const system_type&) const noexcept override {return 0;}

    void initialize(const system_type&) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("connection kind = ", this->kind_,
                         ", number of bonds = ", potentials_.size());
        return;
    }

    void update       (const system_type&)                  override {return;}
    void update_margin(const real_type, const system_type&) override {return;}

    std::string name() const override {return "Dummy:"_s + this->kind_;}

    void write_topology(topology_type& topol) const override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(this->kind_.empty() || this->kind_ == "none")
        {
            MJOLNIR_LOG_WARN("DummyInteraction does not have `topology` field. "
                             "It does not do anything.");
            return;
        }

        for(const auto& idxs : this->potentials_)
        {
            const auto i = idxs[0];
            const auto j = idxs[1];
            topol.add_connection(i, j, this->kind_);
        }
        return;
    }

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

  private:
    connection_kind_type kind_;
    container_type       potentials_;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>

namespace mjolnir
{

// harmonic
extern template class DummyInteraction<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class DummyInteraction<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class DummyInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class DummyInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif// MJOLNIR_INTERACTION_DUMMY_INTERACTION_HPP
