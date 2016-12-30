#ifndef MJOLNIR_EXCLUDED_VOLUME_POTENTIAL
#define MJOLNIR_EXCLUDED_VOLUME_POTENTIAL
#include "GlobalPotentialBase.hpp"
#include <cmath>

namespace mjolnir
{

template<typename traitsT>
class ExcludedVolumePotential : public GlobalPotentialBase<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef real_type parameter_type;

    // rc = 2.5 * sigma
    constexpr static real_type cutoff_ratio = 2.5;

  public:
    ExcludedVolumePotential() = default;
    ~ExcludedVolumePotential() override = default;

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
    
    real_type max_cutoff_length() const;

  private:

    real_type epsilon_;
    std::vector<parameter_type> radii_;
};

template<typename traitsT>
inline void
ExcludedVolumePotential<traitsT>::emplace(const parameter_type& radius)
{
    radii_.push_back(radius);
    return ;
}

template<typename traitsT>
inline void
ExcludedVolumePotential<traitsT>::emplace(parameter_type&& radius)
{
    radii_.emplace_back(std::forward<parameter_type>(radius));
    return ;
}

template<typename traitsT>
inline void
ExcludedVolumePotential<traitsT>::set_radii(const std::vector<parameter_type>& radii)
{
    radii_ = radii;
    return ;
}

template<typename traitsT>
inline void
ExcludedVolumePotential<traitsT>::set_radii(std::vector<parameter_type>&& radii)
{
    radii_ = std::forward<std::vector<parameter_type>>(radii);
    return ;
}


template<typename T>
inline typename T::real_type ExcludedVolumePotential<T>::potential(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    const real_type d = radii_[i] + radii_[j];
    const real_type d_r = d / r;
    const real_type dr3 = d_r * d_r * d_r;
    const real_type dr6 = dr3 * dr3;
    const real_type dr12 = dr6 * dr6;
    return epsilon_ * dr12;
}
 
template<typename T>
inline typename T::real_type ExcludedVolumePotential<T>::derivative(
        const std::size_t i, const std::size_t j, const real_type r) const
{
    const real_type d = radii_[i] + radii_[j];
    const real_type rinv = 1. / r;
    const real_type d_r = d * rinv;
    const real_type dr3 = d_r * d_r * d_r;
    const real_type dr6 = dr3 * dr3;
    const real_type dr12 = dr6 * dr6;
    return -12 * epsilon_ * dr12 * rinv;
}

} // mjolnir
#endif /* MJOLNIR_EXCLUDED_VOLUME_POTENTIAL */
