#ifndef JARNGREIPR_GEOMETRY_ANGLE
#define JARNGREIPR_GEOMETRY_ANGLE
#include "distance.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <mjolnir/util/scalar_type_of.hpp>

namespace jarngreipr
{

template<typename coordT>
inline typename scalar_type_of<coordT>::type
cos_theta(const coordT& lhs, const coordT& rhs)
{
    const auto dist2     = mjolnir::length_sq(lhs) * mjolnir::length_sq(rhs);
    const auto dot       = mjolnir::dot_product(lhs, rhs);
    const auto cos_theta = dot * mjolnir::fast_inv_sqrt(dist2);
         if(cos_theta < -1.) return -1.;
    else if(cos_theta >  1.) return 1.;
    else                     return cos_theta;
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
cos_theta(const coordT& p1, const coordT& p2, const coordT& p3)
{
    return cos_theta(p1 - p2, p3 - p2);
}

template<typename coordT>
typename scalar_type_of<coordT>::type
angle(const coordT& lhs, const coordT& rhs)
{
    return std::acos(cos_theta(lhs, rhs));
}

template<typename traitsT>
typename traitsT::real_type
angle(const coordT& p1, const coordT& p2, const coordT& p3)
{
    return std::acos(cos_theta(p1, p2, p3));
}

} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_ANGLE */
