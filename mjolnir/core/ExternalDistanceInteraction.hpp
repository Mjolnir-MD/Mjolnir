#ifndef MJOLNIR_BOX_INTEARACTION_BASE
#define MJOLNIR_BOX_INTEARACTION_BASE
#include <mjolnir/core/ExternalInteractionBase.hpp>

namespace mjolnir
{

/*! @brief Interaction between particle and a region based on their distance. *
 *  @details shapeT represents the shape. It provides a method to calculate   *
 *           distance between particle and the shape, force direction, and    *
 *           neighbor-list.                                                   */
template<typename traitsT, typename potentialT, typename shapeT>
class ExternalDistanceInteraction final : public ExternalForceInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef shapeT     shape_type;
    typedef ExternalInteractionBase<traits_type> base_type;
    typedef typename base_type::real_type        real_type;
    typedef typename base_type::coordinate_type  coordinate_type;
    typedef typename base_type::system_type      system_type;
    typedef typename base_type::boundary_type    boundary_type;
    typedef typename base_type::particle_type    particle_type;

  public:

    ~ExternalDistanceInteraction() override = default;

    ExternalDistanceInteraction(shape_type&& shape, potential_type&& pot)
        : shape_(std::move(shape)), potential_(std::move(pot))
    {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const real_type dt) override
    {
        this->shape_.set_cutoff(potential_.max_cutoff_length());
        this->shape_.initialize(sys);
        this->shape_.update(sys);
    }

    /*! @brief update parameters (e.g. dt, temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing),     *
     *           the method is bound to call this function after changing     *
     *           parameters.                                                  */
    void reconstruct(const system_type& sys, const real_type dt) override
    {
        this->potential_.update(sys);
        this->shape_.set_cutoff(potential_.max_cutoff_length());
        this->shape_.update(sys);
    }

    //! @brief calculate force, update spatial partition (reduce mergin) inside.
    void      calc_force (system_type&)             override;
    //! @brief calculate energy, do nothing else.
    real_type calc_energy(system_type const&) const override;

    std::string name() const noexcept {return "ExternalDistance";}

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
        const real_type dist =
            this->shape_.calc_distance(sys[i].position, sys.boundary());

        const real_type dV = potential_.derivative(i, dist);
        if(dV == 0.0){continue;}

        sys[i].force +=
            -dV * shape_.calc_force_direction(sys[i].position, sys.boundary());
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
real_type ExternalDistanceInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const
{
    real_type E = 0.0;
    for(std::size_t i : this->shape_.neighbors())
    {
        const real_type d =
            this->shape_.calc_distance(sys[i].position, sys.boundary());

        E += this->potential_.potential(i, d);
    }
    return E;
}

} // mjolnir
#endif//MJOLNIR_BOX_INTEARACTION_BASE
