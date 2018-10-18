#ifndef MJOLNIR_BOND_LENGTH_INTERACTION
#define MJOLNIR_BOND_LENGTH_INTERACTION
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <cmath>
#include <iostream>

namespace mjolnir
{

/*! @brief calculate energy and force of Bond length type local interaction */
template<typename traitsT, typename potentialT>
class BondLengthInteraction : public LocalInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef LocalInteractionBase<traits_type>   base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::topology_type   topology_type;
    typedef typename base_type::connection_kind_type connection_kind_type;

    typedef std::array<std::size_t, 2>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    BondLengthInteraction(const connection_kind_type kind,
                          const container_type& pot)
        : kind_(kind), potentials(pot)
    {}
    BondLengthInteraction(const connection_kind_type kind,
                          container_type&& pot)
        : kind_(kind), potentials(std::move(pot))
    {}
    ~BondLengthInteraction() override = default;

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(BondLengthInteraction::initialize(), 0);
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of bonds = ", potentials.size());
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials)
        {
            item.second.update(sys);
        }
    }

    std::string name() const override
    {return "BondLength:"_s + potential_type::name();}

    void write_topology(topology_type&) const override;

  private:
    connection_kind_type kind_;
    container_type potentials;
};

template<typename traitsT, typename potentialT>
void BondLengthInteraction<traitsT, potentialT>::calc_force(
        system_type& sys) const noexcept
{
    for(const auto& idxp : this->potentials)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];

        const auto dpos =
            sys.adjust_direction(sys[idx1].position - sys[idx0].position);

        const real_type len2 = length_sq(dpos); // l^2
        const real_type rlen = rsqrt(len2);     // 1/l
        const real_type force = -1 * idxp.second.derivative(len2 * rlen);
        // here, L^2 * (1 / L) = L.

        const coordinate_type f = dpos * (force * rlen);
        sys[idx0].force -= f;
        sys[idx1].force += f;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename BondLengthInteraction<traitsT, potentialT>::real_type
BondLengthInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.;
    for(const auto& idxp : this->potentials)
    {
        E += idxp.second.potential(length(sys.adjust_direction(
                sys[idxp.first[1]].position - sys[idxp.first[0]].position)));
    }
    return E;
}

template<typename traitsT, typename potentialT>
void BondLengthInteraction<traitsT, potentialT>::write_topology(
        topology_type& topol) const
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

} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
