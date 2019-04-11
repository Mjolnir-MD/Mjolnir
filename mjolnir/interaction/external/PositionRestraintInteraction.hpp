#ifndef MJOLNIR_INTERACTION_EXTERNAL_POSITION_RESTRAINT_INTERACTION_HPP
#define MJOLNIR_INTERACTION_EXTERNAL_POSITION_RESTRAINT_INTERACTION_HPP
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>

namespace mjolnir
{

//! @brief Interaction between particles and points. 
template<typename traitsT, typename potentialT>
class PositionRestraintInteraction final
    : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using potential_type  = potentialT;
    using base_type       = ExternalForceInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

    // list of {pair of a particle-id and a position at which
    //          the particle would be restrained}
    using shape_type = std::vector<std::pair<std::size_t, coordinate_type>>;

  public:

    PositionRestraintInteraction(
        shape_type shape, potential_type pot)
        : shape_(std::move(shape)), potential_(std::move(pot))
    {}
    ~PositionRestraintInteraction() override = default;

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(system_type const&) const noexcept override;

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        // update system parameters such as temperature, ionic-strength, etc.,
        this->potential_.update(sys);
        return;
    }

    void update(const system_type& sys) override
    {
        // update system parameters such as temperature, ionic-strength, etc.
        this->potential_.update(sys);
        return;
    }

    // no margin to be updated exists.
    void update_margin(const real_type, const system_type&) override {return;}

    std::string name() const override
    {return "PositionRestraint:"_s + potential_.name();}

  private:

    shape_type     shape_;
    potential_type potential_;
};

template<typename traitsT, typename potT>
void PositionRestraintInteraction<traitsT, potT>::calc_force(
        system_type& sys) const noexcept
{
    for(const auto& pidpos : this->shape_)
    {
        const auto  pid = pidpos.first;  // particle index
        const auto& pos = pidpos.second; // restrained position

        const auto& r_i = sys.position(pid);
        const auto  dr  = sys.adjust_direction(pos - r_i);

        const auto dist = math::length(dr);
        const auto dV   = this->potential_.derivative(pid, dist);
        if(dV == 0.0){continue;}

        sys.force(pid) += dV * dr;
    }
    return ;
}

template<typename traitsT, typename potT>
typename PositionRestraintInteraction<traitsT, potT>::real_type
PositionRestraintInteraction<traitsT, potT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const auto& pidpos : this->shape_)
    {
        const auto  pid = pidpos.first;  // particle index
        const auto& pos = pidpos.second; // restrained position

        const auto& r_i = sys.position(pid);
        const auto  dr  = sys.adjust_direction(pos - r_i);

        const auto dist = math::length(dr);
        E += this->potential_.potential(pid, dist);
    }
    return E;
}

} // mjolnir
#endif//MJOLNIR_BOX_INTEARACTION_BASE
