#ifndef JARNGREIPR_GEOMETRY_DISTANCE
#define JARNGREIPR_GEOMETRY_DISTANCE
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <mjolnir/math/Vector.hpp>

namespace jarngreipr
{

template<typename coordT>
inline typename mjolnir::scalar_type_of<coordT>::type
distance_sq(const coordT& lhs, const coordT& rhs)
{
    return mjolnir::length_sq(lhs - rhs);
}

template<typename coordT>
inline typename mjolnir::scalar_type_of<coordT>::type
distance(const coordT& lhs, const coordT& rhs)
{
    return std::sqrt(distance_sq(lhs, rhs));
}


template<typename traitsT>
inline typename traitsT::real_type
distance_sq(const PDBAtom<traitsT>& lhs, const PDBAtom<traitsT>& rhs)
{
    return mjolnir::length_sq(lhs.position - rhs.position);
}

template<typename traitsT>
inline typename traitsT::real_type
distance(const PDBAtom<traitsT>& lhs, const PDBAtom<traitsT>& rhs)
{
    return std::sqrt(distance_sq(lhs, rhs));
}

template<typename traitsT>
typename traitsT::real_type
min_distance_sq(const PDBResidue<traitsT>& lhs, const PDBResidue<traitsT>& rhs)
{
    auto min_dist = std::numeric_limits<typename traitsT::real_type>::max();
    for(auto iter = lhs.cbegin(); iter != lhs.cend(); ++iter)
    {
        for(auto jter = rhs.cbegin(); jter != rhs.cend(); ++jter)
        {
            const auto dist = distance_sq(*iter, *jter);
            if(min_dist > dist) min_dist = dist;
        }
    }
    return min_dist;
}

template<typename traitsT>
inline typename traitsT::real_type
min_distance(const PDBResidue<traitsT>& lhs, const PDBResidue<traitsT>& rhs)
{
    return std::sqrt(min_distance_sq(lhs, rhs));
}

template<typename traitsT>
typename traitsT::real_type
min_distance_sq(const PDBChain<traitsT>& lhs, const PDBChain<traitsT>& rhs)
{
    auto min_dist = std::numeric_limits<typename traitsT::real_type>::max();
    for(auto iter = lhs.cbegin(); iter != lhs.cend(); ++iter)
    {
        for(auto jter = rhs.cbegin(); jter != rhs.cend(); ++jter)
        {
            const auto dist = distance_sq(*iter, *jter);
            if(min_dist > dist) min_dist = dist;
        }
    }
    return min_dist;
}

template<typename traitsT>
inline typename traitsT::real_type
min_distance(const PDBChain<traitsT>& lhs, const PDBChain<traitsT>& rhs)
{
    return std::sqrt(min_distance_sq(lhs, rhs));
}

} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_DISTANCE */
