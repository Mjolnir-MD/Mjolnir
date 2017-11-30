#ifndef JARNGREIPR_GEOMETRY_ANGLE
#define JARNGREIPR_GEOMETRY_ANGLE
#include <mjolnir/math/rsqrt.hpp>
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/util/scalar_type_of.hpp>
#include <cmath>

namespace mjolnir
{

template<typename coordT>
inline typename scalar_type_of<coordT>::type
cos_theta(const coordT& lhs, const coordT& rhs)
{
    const auto dist2     = length_sq(lhs) * length_sq(rhs);
    const auto dot       = dot_product(lhs, rhs);
    const auto cos_theta = dot * mjolnir::rsqrt(dist2);
         if(cos_theta < -1.0){return -1.0;}
    else if(cos_theta >  1.0){return  1.0;}
    else                     {return cos_theta;}
}

template<typename coordT>
inline typename mjolnir::scalar_type_of<coordT>::type
cos_theta(const coordT& p1, const coordT& p2, const coordT& p3)
{
    return cos_theta(p1 - p2, p3 - p2);
}

template<typename coordT>
typename mjolnir::scalar_type_of<coordT>::type
angle(const coordT& lhs, const coordT& rhs)
{
    return std::acos(cos_theta(lhs, rhs));
}

template<typename coordT>
typename mjolnir::scalar_type_of<coordT>::type
angle(const coordT& p1, const coordT& p2, const coordT& p3)
{
    return std::acos(cos_theta(p1, p2, p3));
}

template<typename coordT>
typename mjolnir::scalar_type_of<coordT>::type
angle(const std::unique_ptr<Bead<coordT>>& lhs,
      const std::unique_ptr<Bead<coordT>>& rhs)
{
    return std::acos(cos_theta(lhs->position(), rhs->position()));
}

template<typename coordT>
typename mjolnir::scalar_type_of<coordT>::type
angle(const std::unique_ptr<Bead<coordT>>& p1,
      const std::unique_ptr<Bead<coordT>>& p2,
      const std::unique_ptr<Bead<coordT>>& p3)
{
    return std::acos(cos_theta(p1->position(), p2->position(), p3->position()));
}

} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_ANGLE */
