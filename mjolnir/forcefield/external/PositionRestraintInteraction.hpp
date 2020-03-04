#ifndef MJOLNIR_INTERACTION_EXTERNAL_POSITION_RESTRAINT_INTERACTION_HPP
#define MJOLNIR_INTERACTION_EXTERNAL_POSITION_RESTRAINT_INTERACTION_HPP
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <tuple>

namespace mjolnir
{

// Interaction between particles and fixed points.
//
// XXX
// It re-uses local potentials that can be found in mjolnir/potential/local/*,
// e.g. HarmonicPotential, GaussianPotential, and so on.
//
template<typename traitsT, typename local_potentialT>
class PositionRestraintInteraction final
    : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using potential_type  = local_potentialT;
    using base_type       = ExternalForceInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

    using potential_container_type = std::vector<
            //        {particle-id, fixed point,    potential class}
            std::tuple<std::size_t, coordinate_type, potential_type>
        >;

  public:

    explicit PositionRestraintInteraction(potential_container_type&& pots)
        : potentials_(std::move(pots))
    {}
    ~PositionRestraintInteraction() override {}

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(system_type const&) const noexcept override;

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        // update system parameters such as temperature, ionic-strength, etc.,
        for(auto& pot : potentials_)
        {
            std::get<2>(pot).update(sys);
        }
        return;
    }

    void update(const system_type& sys) override
    {
        // update system parameters such as temperature, ionic-strength, etc.
        for(auto& pot : potentials_)
        {
            std::get<2>(pot).update(sys);
        }
        return;
    }

    // no margin to be updated exists.
    void reduce_margin(const real_type, const system_type&) override {return;}
    void  scale_margin(const real_type, const system_type&) override {return;}

    std::string name() const override
    {return "PositionRestraint:"_s + potential_type::name();}

    // for tests
    potential_container_type const& potentials() const noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new PositionRestraintInteraction(
                potential_container_type(potentials_));
    }

  private:

    potential_container_type potentials_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potT>
void PositionRestraintInteraction<traitsT, potT>::calc_force(
        system_type& sys) const noexcept
{
    for(const auto& pots : this->potentials_)
    {
        const auto  pid = std::get<0>(pots); // particle index
        const auto& pos = std::get<1>(pots); // restrained position
        const auto& pot = std::get<2>(pots); // potential form

        const auto& r_i = sys.position(pid);
        const auto  dr  = sys.adjust_direction(pos - r_i);

        const auto rlen = math::rlength(dr); // 1 / |dr|
        const auto dist = math::length_sq(dr) * rlen;
        const auto dV   = pot.derivative(dist);
        if(dV == 0.0){continue;}

        sys.force(pid) += (dV * rlen) * dr;
    }
    return ;
}

template<typename traitsT, typename potT>
typename PositionRestraintInteraction<traitsT, potT>::real_type
PositionRestraintInteraction<traitsT, potT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const auto& pots : this->potentials_)
    {
        const auto  pid = std::get<0>(pots); // particle index
        const auto& pos = std::get<1>(pots); // restrained position
        const auto& pot = std::get<2>(pots); // potential form

        const auto& r_i = sys.position(pid);
        const auto  dr  = sys.adjust_direction(pos - r_i);

        const auto dist = math::length(dr);
        E += pot.potential(dist);
    }
    return E;
}

} // mjolnir
#endif//MJOLNIR_BOX_INTEARACTION_BASE
