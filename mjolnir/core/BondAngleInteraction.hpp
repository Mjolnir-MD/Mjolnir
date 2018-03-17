#ifndef MJOLNIR_CORE_BOND_ANGLE_INTERACTION
#define MJOLNIR_CORE_BOND_ANGLE_INTERACTION
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/constants.hpp>
#include <mjolnir/math/rsqrt.hpp>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename potentialT>
class BondAngleInteraction : public LocalInteractionBase<traitsT>
{
  public:
    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef LocalInteractionBase<traits_type>   base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::particle_type   particle_type;
    typedef std::array<std::size_t, 3>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    BondAngleInteraction() = default;
    ~BondAngleInteraction() = default;

    BondAngleInteraction(const container_type& pot): potentials(pot){}
    BondAngleInteraction(container_type&& pot): potentials(std::move(pot)){}

    void      calc_force (system_type&)        const noexcept override;
    real_type calc_energy(const system_type& ) const noexcept override;

    void update(const system_type& sys, const real_type dt) override
    {
        for(auto& item : potentials)
        {
            item.second.update(sys, dt);
        }
    }

    std::size_t size() const noexcept {return potentials.size();}

    const_iterator begin()  const noexcept {return potentials.begin();}
    const_iterator end()    const noexcept {return potentials.end();}
    const_iterator cbegin() const noexcept {return potentials.begin();}
    const_iterator cend()   const noexcept {return potentials.end();}

    void append(std::unique_ptr<LocalInteractionBase<traitsT>>&& other) override;

    std::string name() const override
    {return "BondAngle:" + potentials.front().second.name();}

  private:
    container_type potentials;
};

template<typename traitsT, typename potentialT>
void
BondAngleInteraction<traitsT, potentialT>::calc_force(system_type& sys) const noexcept
{
    for(const auto& idxp : this->potentials)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];
        const std::size_t idx2 = idxp.first[2];

        const coordinate_type r_ij =
            sys.adjust_direction(sys[idx0].position - sys[idx1].position);

        const real_type       inv_len_r_ij = rsqrt(length_sq(r_ij));
        const coordinate_type r_ij_reg     = r_ij * inv_len_r_ij;

        const coordinate_type r_kj =
            sys.adjust_direction(sys[idx2].position - sys[idx1].position);

        const real_type       inv_len_r_kj = rsqrt(length_sq(r_kj));
        const coordinate_type r_kj_reg     = r_kj * inv_len_r_kj;

        const real_type dot_ijk = dot_product(r_ij_reg, r_kj_reg);
        const real_type cos_theta = (-1. <= dot_ijk && dot_ijk <= 1.)
            ? dot_ijk : std::copysign(1.0, dot_ijk); // for numerical robustness

        const real_type theta = std::acos(cos_theta);
        const real_type coef  = -(idxp.second.derivative(theta));

        const real_type sin_theta    = std::sin(theta);
        const real_type coef_inv_sin =
            (sin_theta > constants<real_type>::tolerance) ?
            coef / sin_theta : coef / constants<real_type>::tolerance;

        const coordinate_type Fi =
            (coef_inv_sin * inv_len_r_ij) * (cos_theta * r_ij_reg - r_kj_reg);

        const coordinate_type Fk =
            (coef_inv_sin * inv_len_r_kj) * (cos_theta * r_kj_reg - r_ij_reg);

        sys[idx0].force += Fi;
        sys[idx1].force -= (Fi + Fk);
        sys[idx2].force += Fk;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename BondAngleInteraction<traitsT, potentialT>::real_type
BondAngleInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.;
    for(const auto& idxp : this->potentials)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];
        const std::size_t idx2 = idxp.first[2];

        const coordinate_type v_2to1 =
            sys.adjust_direction(sys[idx0].position - sys[idx1].position);
        const coordinate_type v_2to3 =
            sys.adjust_direction(sys[idx2].position - sys[idx1].position);

        const real_type lensq_v21   = length_sq(v_2to1);
        const real_type lensq_v23   = length_sq(v_2to3);
        const real_type dot_v21_v23 = dot_product(v_2to1, v_2to3);

        const real_type dot_ijk   = dot_v21_v23 * rsqrt(lensq_v21 * lensq_v23);
        const real_type cos_theta = (-1. <= dot_ijk && dot_ijk <= 1.)
                                    ? dot_ijk : std::copysign(1.0, dot_ijk);
        const real_type theta     = std::acos(cos_theta);

        E += idxp.second.potential(theta);
    }
    return E;
}

template<typename traitsT, typename potentialT>
void BondAngleInteraction<traitsT, potentialT>::append(
        std::unique_ptr<LocalInteractionBase<traitsT>>&& other)
{
    const BondAngleInteraction<traitsT, potentialT>* rptr =
        dynamic_cast<BondAngleInteraction<traitsT, potentialT>*>(other.get());
    if(rptr == nullptr)
    {
        throw std::invalid_argument("mjolnir::BondAngleInteraction::append: "
                "non-subclass appears!");
    }
    this->potentials.reserve(this->potentials.size() + rptr->size());
    std::copy(rptr->begin(), rptr->end(), std::back_inserter(this->potentials));
    other.release();
    return;
}

}// mjolnir
#endif /* MJOLNIR_BOND_ANGL_INTERACTION */
