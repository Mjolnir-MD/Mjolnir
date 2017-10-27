#ifndef JARNGREIPR_GEOMETRY_DISTANCE
#define JARNGREIPR_GEOMETRY_DISTANCE
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <mjolnir/math/Vector.hpp>

namespace mjolnir
{

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance_sq(const Bead<coordT>& lhs, const Bead<coordT>& rhs) noexcept
{
    return distance_sq(lhs.position(), rhs.position());
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance(const Bead<coordT>& lhs, const Bead<coordT>& rhs) noexcept
{
    return std::sqrt(distance_sq(lhs, rhs));
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance_sq(const PDBAtom<coordT>& lhs, const PDBAtom<coordT>& rhs) noexcept
{
    return distance_sq(lhs.position, rhs.position);
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance(const PDBAtom<coordT>& lhs, const PDBAtom<coordT>& rhs) noexcept
{
    return std::sqrt(distance_sq(lhs, rhs));
}

/* ------------------------------ min_distance ------------------------------ */

// TODO: refactoring (according to DRY)

template<typename coordT>
typename scalar_type_of<coordT>::type
min_distance_sq(const std::vector<PDBAtom<coordT>>& lhs,
                const std::vector<PDBAtom<coordT>>& rhs)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    auto min_dist = std::numeric_limits<real_type>::infinity();

    for(auto l : lhs)
    {
        for(auto r : rhs)
        {
            const auto dist = distance_sq(l, r);
            if(min_dist > dist) {min_dist = dist;}
        }
    }
    return min_dist;
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
min_distance(const std::vector<PDBAtom<coordT>>& lhs,
             const std::vector<PDBAtom<coordT>>& rhs)
{
    return std::sqrt(min_distance_sq(lhs, rhs));
}

template<typename coordT>
typename scalar_type_of<coordT>::type
min_distance_sq(const PDBResidue<coordT>& lhs, const PDBResidue<coordT>& rhs)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    auto min_dist = std::numeric_limits<real_type>::infinity();

    for(auto l : lhs)
    {
        for(auto r : rhs)
        {
            const auto dist = distance_sq(l, r);
            if(min_dist > dist) {min_dist = dist;}
        }
    }
    return min_dist;
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
min_distance(const PDBResidue<coordT>& lhs, const PDBResidue<coordT>& rhs)
{
    return std::sqrt(min_distance_sq(lhs, rhs));
}

template<typename coordT>
typename scalar_type_of<coordT>::type
min_distance_sq(const PDBChain<coordT>& lhs, const PDBChain<coordT>& rhs)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    auto min_dist = std::numeric_limits<real_type>::infinity();

    for(auto l : lhs)
    {
        for(auto r : rhs)
        {
            const auto dist = distance_sq(l, r);
            if(min_dist > dist) {min_dist = dist;}
        }
    }
    return min_dist;
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
min_distance(const PDBChain<coordT>& lhs, const PDBChain<coordT>& rhs)
{
    return std::sqrt(min_distance_sq(lhs, rhs));
}

/* ---------------------------- distance_center ----------------------------- */
// TODO: refactoring (according to DRY)

template<typename coordT>
typename scalar_type_of<coordT>::type
distance_center_sq(const std::vector<PDBAtom<coordT>>& lhs,
                   const std::vector<PDBAtom<coordT>>& rhs) noexcept
{
    return distance_sq(center(lhs), center(rhs));
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance_center(const std::vector<PDBAtom<coordT>>& lhs,
                const std::vector<PDBAtom<coordT>>& rhs) noexcept
{
    return std::sqrt(distance_center_sq(lhs, rhs));
}

template<typename coordT>
typename scalar_type_of<coordT>::type
distance_center_sq(const PDBResidue<coordT>& lhs,
                   const PDBResidue<coordT>& rhs) noexcept
{
    return distance_sq(center(lhs), center(rhs));
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance_center(const PDBResidue<coordT>& lhs,
                const PDBResidue<coordT>& rhs) noexcept
{
    return std::sqrt(distance_center_sq(lhs, rhs));
}

template<typename coordT>
typename scalar_type_of<coordT>::type
distance_center_sq(const PDBChain<coordT>& lhs,
                   const PDBChain<coordT>& rhs) noexcept
{
    return distance_sq(center(lhs), center(rhs));
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance_center(const PDBChain<coordT>& lhs,
                const PDBChain<coordT>& rhs) noexcept
{
    return std::sqrt(distance_center_sq(lhs, rhs));
}

} // jarngreipr
#endif /* JARNGREIPR_GEOMETRY_DISTANCE */
