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
template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class BondLengthInteraction : public LocalInteractionBase<traitsT, 2>
{
  public:

    typedef traitsT traits_type;
    typedef LocalInteractionBase<traits_type, 2> base_type;
    typedef typename base_type::time_type        time_type;
    typedef typename base_type::real_type        real_type;
    typedef typename base_type::coordinate_type  coordinate_type;
    typedef typename base_type::particle_type    particle_type;
    typedef typename base_type::particle_ptrs    particle_ptrs;
    typedef LocalPotentialBase<traits_type>      potential_type;

  public:

    BondLengthInteraction() = default;
    BondLengthInteraction(const boundaryT& b): boundary(b){};
    ~BondLengthInteraction() override = default;

    void
    calc_force(particle_ptrs ps,
               const potential_type& pot) const override;

    real_type
    calc_energy(const particle_ptrs ps,
                const potential_type& pot) const override;

  private:

    boundaryT boundary;
};

template<typename traitsT, typename boundaryT>
void BondLengthInteraction<traitsT, boundaryT>::calc_force(
        particle_ptrs ps, const potential_type& pot) const
{
    const coordinate_type dpos = boundary(ps[1]->position - ps[0]->position);
    const real_type len = length(dpos);
    const real_type f = -1 * pot.derivative(len);

    if(std::abs(f) < std::numeric_limits<real_type>::epsilon()) return;

    const coordinate_type force = dpos * (f / len);
    ps[0]->force -= force;
    ps[1]->force += force;
    return;
}

template<typename traitsT, typename boundaryT>
typename BondLengthInteraction<traitsT, boundaryT>::real_type
BondLengthInteraction<traitsT, boundaryT>::calc_energy(
        const particle_ptrs ps, const potential_type& pot) const
{
    return pot.potential(length(boundary(ps[0]->position - ps[1]->position)));
}

} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
