#ifndef MJOLNIR_BOND_LENGTH_INTERACTION
#define MJOLNIR_BOND_LENGTH_INTERACTION
#include "LocalInteractionBase.hpp"
#include "BoundaryCondition.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <limits>
#include <cmath>

namespace mjolnir
{

/*! @brief calculate energy and force of Bond length type local interaction */
template<typename traitsT, typename potentialT,
         typename boundaryT = UnlimitedBoundary<traitsT>>
class BondLengthInteraction : public LocalInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef boundaryT  boundary_type;
    typedef LocalInteractionBase<traits_type>   base_type;
    typedef typename base_type::time_type       time_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::particle_type   particle_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef std::array<std::size_t, 2>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;

  public:

    BondLengthInteraction() = default;
    ~BondLengthInteraction() override = default;

    BondLengthInteraction(const container_type& pot): potentials(pot){}
    BondLengthInteraction(container_type&& pot)
        : potentials(std::forward<container_type>(pot)){}
    BondLengthInteraction(const BondLengthInteraction&) = default;
    BondLengthInteraction(BondLengthInteraction&&)      = default;
    BondLengthInteraction& operator=(const BondLengthInteraction&) = default;
    BondLengthInteraction& operator=(BondLengthInteraction&&)      = default;

    void
    calc_force(particle_container_type& pcon) const override;

    real_type
    calc_energy(const particle_container_type& pcon) const override;

    void
    reset_parameter(const std::string& name, const real_type val) override;

  private:

    void
    calc_force(particle_type& p1, particle_type& p2,
               const potential_type& pot) const;
    real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const potential_type& pot) const;
  private:
    container_type potentials;
};

template<typename traitsT, typename potentialT, typename boundaryT>
void BondLengthInteraction<traitsT, potentialT, boundaryT>::calc_force(
        particle_container_type& pcon) const
{
    for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
        this->calc_force(pcon[iter->first[0]], pcon[iter->first[1]],
                         iter->second);
    return;
}

template<typename traitsT, typename potentialT, typename boundaryT>
typename BondLengthInteraction<traitsT, potentialT, boundaryT>::real_type
BondLengthInteraction<traitsT, potentialT, boundaryT>::calc_energy(
        const particle_container_type& pcon) const
{
    real_type E = 0.;
    for(auto iter = potentials.cbegin(); iter != potentials.cend(); ++iter)
        E += this->calc_energy(pcon[iter->first[0]], pcon[iter->first[1]],
                               iter->second);
    return E;
}

template<typename traitsT, typename potentialT, typename boundaryT>
void BondLengthInteraction<traitsT, potentialT, boundaryT>::calc_force(
        particle_type& p1, particle_type& p2, const potential_type& pot) const
{
    const coordinate_type dpos =
        boundary_type::adjust_direction(p2.position - p1.position);
    const real_type len = length(dpos);
    const real_type f = -1 * pot.derivative(len);

    if(std::abs(f) < std::numeric_limits<real_type>::epsilon()) return;

    const coordinate_type force = dpos * (f / len);
    p1.force -= force;
    p2.force += force;
    return;
}

template<typename traitsT, typename potentialT, typename boundaryT>
typename BondLengthInteraction<traitsT, potentialT, boundaryT>::real_type
BondLengthInteraction<traitsT, potentialT, boundaryT>::calc_energy(
        const particle_type& p1, const particle_type& p2,
        const potential_type& pot) const
{
    return pot.potential(length(
                boundary_type::adjust_direction(p1.position - p2.position)));
}


template<typename traitsT, typename potentialT, typename boundaryT>
void BondLengthInteraction<traitsT, potentialT, boundaryT>::reset_parameter(
        const std::string& name, const real_type val)
{
    for(auto iter = potentials.begin(); iter != potentials.end(); ++iter)
        iter->second.reset_parameter(name, val);
    return;
}

} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
