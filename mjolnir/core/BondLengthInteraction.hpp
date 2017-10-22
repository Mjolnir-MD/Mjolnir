#ifndef MJOLNIR_BOND_LENGTH_INTERACTION
#define MJOLNIR_BOND_LENGTH_INTERACTION
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/math/rsqrt.hpp>
#include <cmath>

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
    typedef typename base_type::particle_type   particle_type;
    typedef std::array<std::size_t, 2>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;

  public:

    BondLengthInteraction() = default;
    ~BondLengthInteraction() override = default;

    BondLengthInteraction(const container_type& pot): potentials(pot){}
    BondLengthInteraction(container_type&& pot): potentials(std::move(pot)){}

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

  private:
    container_type potentials;
};

template<typename traitsT, typename potentialT>
void
BondLengthInteraction<traitsT, potentialT>::calc_force(system_type& sys) const noexcept
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

} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
