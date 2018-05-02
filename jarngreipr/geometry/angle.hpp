#ifndef JARNGREIPR_GEOMETRY_ANGLE
#define JARNGREIPR_GEOMETRY_ANGLE
#include <mjolnir/math/rsqrt.hpp>
#include <mjolnir/math/Vector.hpp>
#include <cmath>

namespace jarngreipr
{

template<typename realT>
inline realT cos_theta(const mjolnir::Vector<realT, 3>& lhs,
                       const mjolnir::Vector<realT, 3>& rhs)
{
    const auto dist2 = length_sq(lhs) * length_sq(rhs);
    const auto dot   = mjolnir::dot_product(lhs, rhs);
    const auto cos_t = dot * mjolnir::rsqrt(dist2);
    if(cos_t < -1.0) {return -1.0;} else if(cos_t >  1.0) {return 1.0;}
    return cos_t;
}

template<typename realT>
inline realT cos_theta(const mjolnir::Vector<realT, 3>& p1,
                       const mjolnir::Vector<realT, 3>& p2,
                       const mjolnir::Vector<realT, 3>& p3)
{
    /*      p1     *
     *     /       *
     *    /) theta *
     * p2 ----p3   */
    return cos_theta(p1 - p2, p3 - p2);
}

template<typename realT>
inline realT angle(const mjolnir::Vector<realT, 3>& lhs,
                   const mjolnir::Vector<realT, 3>& rhs)
{
    return std::acos(cos_theta(lhs, rhs));
}

template<typename realT>
inline realT angle(const mjolnir::Vector<realT, 3>& p1,
                   const mjolnir::Vector<realT, 3>& p2,
                   const mjolnir::Vector<realT, 3>& p3)
{
    return std::acos(cos_theta(p1, p2, p3));
}

} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_ANGLE */
