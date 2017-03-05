#ifndef MJOLNIR_GLOBAL_EXCLUDED_VOLUME_INTEARACTION
#define MJOLNIR_GLOBAL_EXCLUDED_VOLUME_INTEARACTION
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/potential/ExcludedVolumePotential.hpp>

namespace mjolnir
{

// specialization for excluded volume potential
template<typename traitsT, typename partitionT, typename boundaryT>
class GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
                                partitionT, boundaryT>
    : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT    traits_type;
    typedef ExcludedVolumePotential<traitsT> potential_type;
    typedef partitionT spatial_partition_type;
    typedef boundaryT  boundary_type;
    typedef GlobalInteractionBase<traitsT> base_type;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::particle_container_type particle_container_type;

  public:
    GlobalDistanceInteraction()  = default;
    ~GlobalDistanceInteraction() = default;

    GlobalDistanceInteraction(potential_type&& pot,
                              spatial_partition_type&& space)
        : potential_(std::forward<potential_type>(pot)),
          spatial_partition_(std::forward<spatial_partition_type>(space))
    {}

    void
    initialize(const particle_container_type& pcon, const time_type dt) override;

    void
    calc_force(particle_container_type& pcon) override;

    real_type
    calc_energy(const particle_container_type& pcon) const override;

    void
    reset_parameter(const std::string& name, const real_type val) override;

  private:
    real_type              cutoff2_;
    potential_type         potential_;
    spatial_partition_type spatial_partition_;
};

template<typename traitsT, typename spaceT, typename boundaryT>
void GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
     spaceT, boundaryT>::calc_force(particle_container_type& pcon)
{
    spatial_partition_.update(pcon);
    const real_type coef = 12 * potential_.epsilon();
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spatial_partition_type::index_list const& partners =
            spatial_partition_.partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = boundary_type::adjust_direction(
                    pcon[j].position - pcon[i].position);
            const real_type r2 = length_sq(rij);
            if(cutoff2_ < r2) continue;

            const real_type sgm  = potential_[i] + potential_[j];
            const real_type invr = fast_inv_sqrt(r2);
            const real_type sr   = sgm * invr;
            const real_type sr2  = sr  * sr;
            const real_type sr6  = sr2 * sr2 * sr2;
            const real_type sr12 = sr6 * sr6;

            const coordinate_type f = (coef * sr12 * invr * invr) * rij;
            pcon[i].force += f;
            pcon[j].force -= f;
        }
    }
    return ;
}

template<typename traitsT, typename spaceT, typename boundaryT>
typename GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
         spaceT, boundaryT>::real_type
GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>, spaceT,
    boundaryT>::calc_energy(const particle_container_type& pcon) const
{
    real_type e = 0.0;
    const real_type eps = potential_.epsilon();
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spatial_partition_type::index_list const& partners =
            spatial_partition_.partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = boundary_type::adjust_direction(
                    pcon[j].position - pcon[i].position);
            const real_type sgm  = potential_[i] + potential_[j];

            const real_type r2   = length_sq(rij);
            const real_type sr   = sgm * fast_inv_sqrt(r2);
            const real_type sr3  = sr * sr * sr;
            const real_type sr12 = sr3 * sr3 * sr3 * sr3;

            e += eps * sr12;
        }
    }
    return e;
}

template<typename traitsT, typename spaceT, typename boundaryT>
void GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
     spaceT, boundaryT>::initialize(const particle_container_type& pcon,
                                    const time_type dt)
{
    const auto rc = potential_.max_cutoff_length();
    cutoff2_ = rc * rc;
    this->spatial_partition_.update(pcon, dt);
    return ;
}

template<typename traitsT, typename spaceT, typename boundaryT>
void GlobalDistanceInteraction<traitsT, ExcludedVolumePotential<traitsT>,
     spaceT, boundaryT>::reset_parameter(const std::string& name, const real_type val)
{
    this->potential_.reset_parameter(name, val);
    return ;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_EXCLUDED_VOLUME_INTEARACTION */
