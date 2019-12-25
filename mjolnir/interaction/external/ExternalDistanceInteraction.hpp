#ifndef MJOLNIR_INTERACTION_EXTERNAL_DISTANCE_INTERACTION_HPP
#define MJOLNIR_INTERACTION_EXTERNAL_DISTANCE_INTERACTION_HPP
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>

namespace mjolnir
{

/*! @brief Interaction between particle and a region based on their distance. *
 *  @details shapeT represents the shape. It provides a method to calculate   *
 *           distance between particle and the shape, force direction, and    *
 *           neighbor-list.                                                   */
template<typename traitsT, typename potentialT, typename shapeT>
class ExternalDistanceInteraction final
    : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using potential_type  = potentialT;
    using shape_type      = shapeT;
    using base_type       = ExternalForceInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

  public:

    ExternalDistanceInteraction(shape_type&& shape, potential_type&& pot)
        : shape_(std::move(shape)), potential_(std::move(pot))
    {}
    ~ExternalDistanceInteraction() override {}

    // calculate force, update spatial partition (reduce margin) inside.
    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(system_type const&) const noexcept override;

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        this->potential_.update(sys); // update system parameters
        this->shape_.initialize(sys, this->potential_);
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys) override
    {
        this->potential_.update(sys); // update system parameters
        this->shape_.initialize(sys, this->potential_);
    }

    void update_margin(const real_type dmargin, const system_type& sys) override
    {
        this->shape_.reduce_margin(dmargin, sys);
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->shape_.scale_margin(scale, sys);
    }

    std::string name() const override
    {return "ExternalDistance:"_s + potential_.name();}

    base_type* clone() const override
    {
        return new ExternalDistanceInteraction(
                shape_type(shape_), potential_type(potential_));
    }

  private:

    shape_type     shape_;
    potential_type potential_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potT, typename spaceT>
void ExternalDistanceInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys) const noexcept
{
    for(std::size_t i : this->shape_.neighbors())
    {
        const auto& ri = sys.position(i);

        const real_type dist = this->shape_.calc_distance(ri, sys.boundary());
        const real_type dV   = this->potential_.derivative(i, dist);
        if(dV == 0.0){continue;}

        const auto f = shape_.calc_force_direction(ri, sys.boundary());
        sys.force(i) += -dV * f;
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename ExternalDistanceInteraction<traitsT, potT, spaceT>::real_type
ExternalDistanceInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(std::size_t i : this->shape_.neighbors())
    {
        const auto&    ri = sys.position(i);
        const real_type d = this->shape_.calc_distance(ri, sys.boundary());
        E += this->potential_.potential(i, d);
    }
    return E;
}

} // mjolnir
#endif//MJOLNIR_BOX_INTEARACTION_BASE
