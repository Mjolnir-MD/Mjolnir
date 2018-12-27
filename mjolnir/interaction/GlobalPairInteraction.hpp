#ifndef MJOLNIR_GLOBAL_PAIR_INTEARACTION
#define MJOLNIR_GLOBAL_PAIR_INTEARACTION
#include <mjolnir/core/GlobalInteractionBase.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/util/string.hpp>
#include <memory>

namespace mjolnir
{

template<typename traitsT, typename potentialT, typename partitionT>
class GlobalPairInteraction final : public GlobalInteractionBase<traitsT>
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
    GlobalPairInteraction()  = default;
    ~GlobalPairInteraction() = default;

    GlobalPairInteraction(potential_type&& pot, partition_type&& part)
        : potential_(std::move(pot)), partition_(std::move(part))
    {}

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(GlobalPairInteraction::initialize(), 0);
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->partition_.initialize(sys, this->potential_);
        this->partition_.update(sys, this->potential_);
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_SCOPE(GlobalPairInteraction::update(), 0);
        MJOLNIR_LOG_INFO("potential is ", this->name());
        this->potential_.update(sys);
        // potential update may change the cutoff length!
        this->partition_.reconstruct(sys, this->potential_);
        this->partition_.update(sys, this->potential_);
    }

    void      calc_force (system_type&)             override;
    real_type calc_energy(const system_type&) const override;

    std::string name() const override
    {return "Pair:"_s + potential_type::name();}

  private:

    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT, typename potT, typename spaceT>
void GlobalPairInteraction<traitsT, potT, spaceT>::calc_force(
        system_type& sys)
{
    partition_.update(sys, this->potential_);
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j     = ptnr.index;
            const auto& param = ptnr.parameter();

            const coordinate_type rij =
                sys.adjust_direction(sys[j].position - sys[i].position);
            const real_type l = length(rij);
            const real_type f_mag = potential_.derivative(l, param);

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
typename GlobalPairInteraction<traitsT, potT, spaceT>::real_type
GlobalPairInteraction<traitsT, potT, spaceT>::calc_energy(
        const system_type& sys) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        for(const auto& ptnr : this->partition_.partners(i))
        {
            const auto  j     = ptnr.index;
            const auto& param = ptnr.parameter();

            const real_type l = length(
                sys.adjust_direction(sys[j].position - sys[i].position));
            e += potential_.potential(l, param);
        }
    }
    return e;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_PAIR_INTEARACTION */
