#ifndef MJOLNIR_LENNARD_JONES_POTENTIAL
#define MJOLNIR_LENNARD_JONES_POTENTIAL
#include "GlobalPotentialBase.hpp"
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

namespace mjolnir
{

/*! @brief Lennard-Jones type potential & derivative                       *
 * designed for global force field. so it doesn't have its own parameters. *
 * V(r)  =  4. * epsilon * ((r/sigma)^12 - (r/sigma)^6))                   *
 * dV/dr = 24. * epsilon / r * ((r/sigma)^6 - 2 * (r/sigma)^12)            */
template<typename traitsT>
class LennardJonesPotential: public GlobalPotentialBase<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::pair<real_type, real_type> parameter_type;

  public:
    LennardJonesPotential() = default;
    ~LennardJonesPotential() = default;

    void emplace(const parameter_type& radius);
    void emplace(parameter_type&& radius);

    void set_radii(const std::vector<parameter_type>& radii);
    void set_radii(std::vector<parameter_type>&& radii);

    parameter_type&       operator[](const std::size_t i)       {return radii_[i];}
    parameter_type const& operator[](const std::size_t i) const {return radii_[i];}
    parameter_type&       at(const std::size_t i)       {return radii_.at(i);}
    parameter_type const& at(const std::size_t i) const {return radii_.at(i);}

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const override;

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const override;

  private:

    std::vector<parameter_type> radii_;
};

template<typename traitsT>
inline void
LennardJonesPotential<traitsT>::emplace(const parameter_type& radius)
{
    radii_.push_back(radius);
    return ;
}

template<typename traitsT>
inline void
LennardJonesPotential<traitsT>::emplace(parameter_type&& radius)
{
    radii_.emplace_back(std::forward<parameter_type>(radius));
    return ;
}

template<typename traitsT>
inline void
LennardJonesPotential<traitsT>::set_radii(const std::vector<parameter_type>& radii)
{
    radii_ = radii;
    return ;
}

template<typename traitsT>
inline void
LennardJonesPotential<traitsT>::set_radii(std::vector<parameter_type>&& radii)
{
    radii_ = std::forward<std::vector<parameter_type>>(radii);
    return ;
}

template<typename traitsT>
inline typename LennardJonesPotential<traitsT>::real_type
LennardJonesPotential<traitsT>::potential(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    const real_type sigma   = 0.5 * (radii_[i].first + radii_[j].first);
    const real_type epsilon = std::sqrt(radii_[i].second * radii_[j].second);

    const real_type r1s1   = sigma / r;
    const real_type r3s3   = r1s1 * r1s1 * r1s1;
    const real_type r6s6   = r3s3 * r3s3;
    const real_type r12s12 = r6s6 * r6s6;
    return 4. * epsilon * (r12s12 - r6s6);
}

template<typename traitsT>
inline typename LennardJonesPotential<traitsT>::real_type
LennardJonesPotential<traitsT>::derivative(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    const real_type sigma   = 0.5 * (radii_[i].first + radii_[j].first);
    const real_type epsilon = std::sqrt(radii_[i].second * radii_[j].second);

    const real_type r1s1   = sigma / r;
    const real_type r3s3   = r1s1 * r1s1 * r1s1;
    const real_type r6s6   = r3s3 * r3s3;
    const real_type r12s12 = r6s6 * r6s6;
    return 24. * epsilon * (r6s6 - 2 * r12s12) / r;
}

} // mjolnir

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
