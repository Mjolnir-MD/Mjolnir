#ifndef MJOLNIR_BOX_INTEARACTION_BASE
#define MJOLNIR_BOX_INTEARACTION_BASE
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
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

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef shapeT     shape_type;
    typedef ExternalForceInteractionBase<traits_type> base_type;
    typedef typename base_type::real_type        real_type;
    typedef typename base_type::coordinate_type  coordinate_type;
    typedef typename base_type::system_type      system_type;
    typedef typename base_type::boundary_type    boundary_type;

  public:

    ExternalDistanceInteraction(shape_type&& shape, potential_type&& pot)
        : shape_(std::move(shape)), potential_(std::move(pot))
    {}
    ~ExternalDistanceInteraction() override = default;

    // calculate force, update spatial partition (reduce mergin) inside.
    void      calc_force (system_type&)             override;
    real_type calc_energy(system_type const&) const override;

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const real_type dt) override
    {
        this->potential_.update(sys); // update system parameters
        this->shape_.initialize(sys, this->potential_);
        this->shape_.update(sys);
    }

    /*! @brief update parameters (e.g. dt, temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing),     *
     *           the method is bound to call this function after changing     *
     *           parameters.                                                  */
    void reconstruct(const system_type& sys, const real_type dt) override
    {
        this->potential_.update(sys); // update system parameters
        this->shape_.reconstruct(sys, this->potential_);
        this->shape_.update(sys);
    }

    std::string name() const override
    {return "ExternalDistance:"_str + potential_.name();}

  private:

    shape_type     shape_;
    potential_type potential_;
};

template<typename traitsT, typename potT, typename spaceT>
void ExternalDistanceInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys)
{
    this->shape_.update(sys); // update neighbor list...

    for(std::size_t i : this->shape_.neighbors())
    {
        const auto& ri = sys[i].position;

        const real_type dist = this->shape_.calc_distance(ri, sys.boundary());
        const real_type dV   = this->potential_.derivative(i, dist);
        if(dV == 0.0){continue;}

        const auto f = shape_.calc_force_direction(ri, sys.boundary());
        sys[i].force += -dV * f;
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename ExternalDistanceInteraction<traitsT, potT, spaceT>::real_type
ExternalDistanceInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const
{
    real_type E = 0.0;
    for(std::size_t i : this->shape_.neighbors())
    {
        const auto&    ri = sys[i].position;
        const real_type d = this->shape_.calc_distance(ri, sys.boundary());
        E += this->potential_.potential(i, d);
    }
    return E;
}

} // mjolnir
#endif//MJOLNIR_BOX_INTEARACTION_BASE
