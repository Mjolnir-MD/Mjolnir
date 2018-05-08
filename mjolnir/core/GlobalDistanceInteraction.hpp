#ifndef MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#define MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#include "GlobalInteractionBase.hpp"
#include <memory>

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename partitionT>
class GlobalDistanceInteraction final : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef partitionT partition_type;
    typedef GlobalInteractionBase<traitsT> base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::boundary_type   boundary_type;

  public:
    GlobalDistanceInteraction()  = default;
    ~GlobalDistanceInteraction() = default;

    GlobalDistanceInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys, const real_type dt) override
    {
        this->partition_.initialize(sys, this->potential_);
        this->partition_.update(sys);
    }

    /*! @brief update parameters (e.g. dt, temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing),     *
     *           the method is bound to call this function after changing     *
     *           parameters.                                                  */
    void update(const system_type& sys, const real_type dt) override
    {
        this->potential_.update(sys);
        // potential update may change the cutoff length!
        this->partition_.reconstruct(sys, this->potential_);
    }

    void      calc_force (system_type&)             override;
    real_type calc_energy(const system_type&) const override;

    std::string name() const override
    {return "Distance:"_str + potential_.name();}

  private:

    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT, typename potT, typename spaceT>
void GlobalDistanceInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys)
{
    partition_.update(sys); // reduce margin and reconstruct pair-list if needed
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        for(auto j : this->partition_.partners(i))
        {
            const coordinate_type rij =
                sys.adjust_direction(sys[j].position - sys[i].position);
            const real_type l = length(rij);
            const real_type f_mag = potential_.derivative(i, j, l);

            // if length exceeds cutoff, potential returns just 0.
            if(f_mag == 0.0){continue;}

            const coordinate_type f = rij * (f_mag / l);
            sys[i].force += f;
            sys[j].force -= f;
        }
    }
    return ;
}

template<typename traitsT, typename potT, typename spaceT>
typename GlobalDistanceInteraction<traitsT, potT, spaceT>::real_type
GlobalDistanceInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        for(auto j : this->partition_.partners(i))
        {
            const real_type l = length(
                sys.adjust_direction(sys[j].position - sys[i].position));
            e += potential_.potential(i, j, l);
        }
    }
    return e;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
