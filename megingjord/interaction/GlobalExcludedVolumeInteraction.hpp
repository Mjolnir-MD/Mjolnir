#ifndef MJOLNIR_GLOBAL_EXCLUDED_VOLUME_INTEARACTION
#define MJOLNIR_GLOBAL_EXCLUDED_VOLUME_INTEARACTION
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/potential/ExcludedVolumePotential.hpp>

// specialization for excluded volume potential
namespace mjolnir
{

template<typename traitsT, typename partitionT>
class GlobalDistanceInteraction<
    traitsT, ExcludedVolumePotential<traitsT>, partitionT>
final : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef partitionT partition_type;
    typedef ExcludedVolumePotential<traitsT>    potential_type;
    typedef GlobalInteractionBase<traitsT>      base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::particle_type   particle_type;
    typedef typename base_type::boundary_type   boundary_type;

  public:
    GlobalDistanceInteraction()  = default;
    ~GlobalDistanceInteraction() = default;

    GlobalDistanceInteraction(potential_type&& pot, partition_type&& space)
        : potential_(std::move(pot)), spatial_partition_(std::move(space))
    {}

    void initialize(const system_type& sys, const real_type dt) override
    {
        this->partition_.initialize(sys);
        this->partition_.update(sys, dt);
    }

    void      calc_force (system_type&)             override;
    real_type calc_energy(const system_type&) const override;

  private:

    potential_type potential_;
    partition_type partition_;
};

template<typename traitsT, typename partitionT>
void GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
     partitionT>::calc_force(system_type& sys)
{
    potential_.update(sys);
    partition_.update(sys);

    const real_type coef = -12 * potential_.epsilon();
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const coordinate_type ri = sys[i].position;
        const real_type si = potential_[i];

        for(auto j : partition_.partners(i))
        {
            const coordinate_type rij =
                sys.adjust_direction(sys[j].position - ri);
            const real_type r2 = length_sq(rij);
            if(cutoff2_ < r2) continue;

            const real_type sgm  = si + potential_[j];
            const real_type invr = fast_inv_sqrt(r2);
            const real_type sr   = sgm * invr;
            if(sr < 1. / potential_type::cutoff_ratio) continue;
            const real_type sr3  = sr  * sr * sr;
            const real_type sr6  = sr3 * sr3;
            const real_type sr12 = sr6 * sr6;

            const coordinate_type f = (coef * sr12 * invr * invr) * rij;
            sys[i].force += f;
            sys[j].force -= f;
        }
    }
    return ;
}

template<typename traitsT, typename partitionT>
typename GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
         partitionT>::real_type
GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
    partitionT>::calc_energy(const system_type& sys) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        for(auto j : this->partition_.partners(i))
        {
            const coordinate_type rij =
                sys.adjust_direction(sys[j].position - sys[i].position);
            const real_type sgm  = potential_[i] + potential_[j];

            const real_type r2   = length_sq(rij);
            const real_type sr   = sgm * fast_inv_sqrt(r2);
            const real_type sr3  = sr  * sr * sr;
            const real_type sr6  = sr3 * sr3;
            const real_type sr12 = sr6 * sr6;
            e += sr12;
        }
    }
    return potential_.epsilon() * e;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_EXCLUDED_VOLUME_INTEARACTION */
