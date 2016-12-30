#ifndef MJOLNIR_LENNARD_JONES_POTENTIAL
#define MJOLNIR_LENNARD_JONES_POTENTIAL
#include "GlobalPotentialBase.hpp"
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <algorithm>
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

    // rc = 2.5 * sigma
    constexpr static real_type cutoff_ratio = 2.5;

  public:
    LennardJonesPotential() = default;
    ~LennardJonesPotential() = default;

    void emplace(const parameter_type& radius);
    void emplace(parameter_type&& radius);

    void set_radii(const std::vector<parameter_type>& radii);
    void set_radii(std::vector<parameter_type>&& radii);

    std::size_t size() const {return radii_.size();}
    void resize(const std::size_t i){radii_.resize(i);}
    void reserve(const std::size_t i){radii_.reserve(i);}
    void clear(){radii_.clear();}
    parameter_type&       operator[](const std::size_t i)       {return radii_[i];}
    parameter_type const& operator[](const std::size_t i) const {return radii_[i];}
    parameter_type&       at(const std::size_t i)       {return radii_.at(i);}
    parameter_type const& at(const std::size_t i) const {return radii_.at(i);}

    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const override;

    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const override;

    real_type max_cutoff_length() const;

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
    if(sigma * cutoff_ratio < r) return 0;

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
    if(sigma * cutoff_ratio < r) return 0;

    const real_type epsilon = std::sqrt(radii_[i].second * radii_[j].second);

    const real_type r1s1   = sigma / r;
    const real_type r3s3   = r1s1 * r1s1 * r1s1;
    const real_type r6s6   = r3s3 * r3s3;
    const real_type r12s12 = r6s6 * r6s6;
    return 24. * epsilon * (r6s6 - 2 * r12s12) / r;
}


template<typename traitsT>
inline typename LennardJonesPotential<traitsT>::real_type
LennardJonesPotential<traitsT>::max_cutoff_length() const
{
    const real_type max_sigma = std::max_element(radii_.cbegin(), radii_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs)
            {return lhs.first < rhs.first;})->first;
    return max_sigma * cutoff_ratio;
}

} // mjolnir

#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
