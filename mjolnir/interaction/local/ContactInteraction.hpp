#ifndef MJOLNIR_INTERACTION_CONTACT_INTERACTION_HPP
#define MJOLNIR_INTERACTION_CONTACT_INTERACTION_HPP
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <cmath>
#include <iostream>

namespace mjolnir
{

// ContactInteraction is basically the same as BondLengthInteraction. The only
// difference is that it has a kind of "Neighbor List" to skip broken contacts.
// BondLengthInteraction considers that all the interactions are kept within
// the cutoff range (if a potential has a cutoff). Contrary, ContactInteraction
// considers any of the interactions occasionally become distant than the cutoff
// range. It becomes faster when a number of the contacts are not formed
// throughout the simulation time.
//
// For the following potentials, there are specialized definitions to speedup
// calculation.
// - GoContactPotential
//   - GoContactInteraction.hpp
// By combining potential calculations, we can omit `sqrt()` that is usually
// used to calculate the distance between particles.
//
template<typename traitsT, typename potentialT>
class ContactInteraction final : public LocalInteractionBase<traitsT>
{
  public:
    using traits_type          = traitsT;
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

    ContactInteraction(const connection_kind_type kind,
                       const container_type&      pot,
                       const real_type            margin = 0.5)
        : kind_(kind), potentials_(pot), margin_(margin)
    {}
    ContactInteraction(const connection_kind_type kind,
                       container_type&&           pot,
                       const real_type            margin = 0.5)
        : kind_(kind), potentials_(std::move(pot)), margin_(margin)
    {}
    ~ContactInteraction() override {}

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of bonds = ", potentials_.size());

        this->cutoff_ = std::max_element(potentials_.begin(), potentials_.end(),
            [](const potential_index_pair& lhs, const potential_index_pair& rhs)
            {
                return lhs.second.cutoff() < rhs.second.cutoff();
            })->second.cutoff();
        this->make_list(sys);
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials_)
        {
            item.second.update(sys);
        }

        this->cutoff_ = std::max_element(potentials_.begin(), potentials_.end(),
            [](const potential_index_pair& lhs, const potential_index_pair& rhs)
            {
                return lhs.second.cutoff() < rhs.second.cutoff();
            })->second.cutoff();
        this->make_list(sys);
        return;
    }

    void update_margin(const real_type dmargin, const system_type& sys) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make_list(sys);
        }
        return;
    }

    std::string name() const override
    {return "Contact:"_s + potential_type::name();}

    void write_topology(topology_type&) const override;

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new ContactInteraction(kind_, container_type(potentials_));
    }

  private:

    void make_list(const system_type& sys)
    {
        this->active_contacts_.clear();
        this->active_contacts_.reserve(potentials_.size());

        // absolute length of margin (this->margin_ is a relative length).
        const real_type abs_margin = this->cutoff_ * this->margin_;

        for(std::size_t i=0; i < this->potentials_.size(); ++i)
        {
            const auto& pot = this->potentials_[i];
            const auto pos0 = sys.position(pot.first[0]);
            const auto pos1 = sys.position(pot.first[1]);
            const auto dpos = sys.adjust_direction(pos1 - pos0);
            const auto len2 = math::length_sq(dpos);

            const auto rc = pot.second.cutoff() + abs_margin;
            if(len2 < rc * rc)
            {
                this->active_contacts_.push_back(i);
            }
        }
        this->current_margin_ = this->cutoff_ * this->margin_;
        return;
    }

  private:
    connection_kind_type kind_;
    container_type potentials_;

    // neighbor list stuff
    real_type cutoff_;
    real_type margin_;
    real_type current_margin_;
    std::vector<std::size_t> active_contacts_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potentialT>
void ContactInteraction<traitsT, potentialT>::calc_force(
        system_type& sys) const noexcept
{
    for(const std::size_t active_contact : active_contacts_)
    {
        const auto& idxp = this->potentials_[active_contact];

        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];

        const auto dpos =
            sys.adjust_direction(sys.position(idx1) - sys.position(idx0));

        const real_type len2 = math::length_sq(dpos); // l^2
        const real_type rlen = math::rsqrt(len2);     // 1/l
        const real_type force = -1 * idxp.second.derivative(len2 * rlen);
        // here, L^2 * (1 / L) = L.

        const coordinate_type f = dpos * (force * rlen);
        sys.force(idx0) -= f;
        sys.force(idx1) += f;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename ContactInteraction<traitsT, potentialT>::real_type
ContactInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const std::size_t active_contact : active_contacts_)
    {
        const auto& idxp = this->potentials_[active_contact];
        E += idxp.second.potential(math::length(sys.adjust_direction(
                sys.position(idxp.first[1]) - sys.position(idxp.first[0]))));
    }
    return E;
}

template<typename traitsT, typename potentialT>
void ContactInteraction<traitsT, potentialT>::write_topology(
        topology_type& topol) const
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

} // mjolnir

#include <mjolnir/interaction/local/GoContactInteraction.hpp>

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/potential/local/GoContactPotential.hpp>
#include <mjolnir/potential/local/GaussianPotential.hpp>

namespace mjolnir
{

// gaussian
extern template class ContactInteraction<SimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
extern template class ContactInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
extern template class ContactInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
extern template class ContactInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
