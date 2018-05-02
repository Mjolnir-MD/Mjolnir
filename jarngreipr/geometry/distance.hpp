#ifndef JARNGREIPR_GEOMETRY_DISTANCE
#define JARNGREIPR_GEOMETRY_DISTANCE
#include <jarngreipr/model/Bead.hpp>
#include <mjolnir/math/Vector.hpp>

namespace jarngreipr
{

template<typename realT>
inline realT distance_sq(const mjolnir::Vector<realT, 3>& lhs,
                         const mjolnir::Vector<realT, 3>& rhs) noexcept
{
    return mjolnir::length_sq(lhs - rhs);
}

template<typename realT>
inline realT distance(const mjolnir::Vector<realT, 3>& lhs,
                      const mjolnir::Vector<realT, 3>& rhs) noexcept
{
    return mjolnir::length(lhs - rhs);
}

} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_DISTANCE */
