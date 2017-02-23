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
class BondLengthInteraction : public LocalInteractionBase<traitsT, 2>
{
  public:

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef boundaryT  boundary_type;
    typedef LocalInteractionBase<traits_type, 2> base_type;
    typedef typename base_type::time_type        time_type;
    typedef typename base_type::real_type        real_type;
    typedef typename base_type::coordinate_type  coordinate_type;
    typedef typename base_type::particle_type    particle_type;

  public:

    BondLengthInteraction() = default;
    ~BondLengthInteraction() override = default;

    BondLengthInteraction(const potential_type& pot): potential(pot){}
    BondLengthInteraction(potential_type&& pot)
        : potential(std::forward<potential_type>(pot)){}
    BondLengthInteraction(const BondLengthInteraction&) = default;
    BondLengthInteraction(BondLengthInteraction&&)      = default;
    BondLengthInteraction& operator=(const BondLengthInteraction&) = default;
    BondLengthInteraction& operator=(BondLengthInteraction&&)      = default;

    void
    calc_force(particle_type& p1, particle_type& p2) const override;

    real_type
    calc_energy(const particle_type& p1, const particle_type& p2) const override;

    void
    reset_parameter(const std::string& name, const real_type val) override
    {
        potential.reset_parameter(name, val);
    }
  private:
    potential_type potential;
};

template<typename traitsT, typename potentialT, typename boundaryT>
void BondLengthInteraction<traitsT, potentialT, boundaryT>::calc_force(
        particle_type& p1, particle_type& p2) const
{
    const coordinate_type dpos =
        boundary_type::adjust_direction(p2.position - p1.position);
    const real_type len = length(dpos);
    const real_type f = -1 * potential.derivative(len);

    if(std::abs(f) < std::numeric_limits<real_type>::epsilon()) return;

    const coordinate_type force = dpos * (f / len);
    p1.force -= force;
    p2.force += force;
    return;
}

template<typename traitsT, typename potentialT, typename boundaryT>
typename BondLengthInteraction<traitsT, potentialT, boundaryT>::real_type
BondLengthInteraction<traitsT, potentialT, boundaryT>::calc_energy(
        const particle_type& p1, const particle_type& p2) const
{
    return potential.potential(length(
                boundary_type::adjust_direction(p1.position - p2.position)));
}

} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
