#ifndef JARNGREIPR_GEOMETRY_DISTANCE
#define JARNGREIPR_GEOMETRY_DISTANCE
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <mjolnir/math/Vector.hpp>

namespace mjolnir
{

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance_sq(const std::unique_ptr<Bead<coordT>>& lhs,
            const std::unique_ptr<Bead<coordT>>& rhs) noexcept
{
    return distance_sq(lhs->position(), rhs->position());
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
distance(const std::unique_ptr<Bead<coordT>>& lhs,
         const std::unique_ptr<Bead<coordT>>& rhs) noexcept
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

// preparation for generics...

template<typename coordT>
inline typename scalar_type_of<coordT>::type
min_distance_sq(const std::unique_ptr<Bead<coordT>>& lhs,
                const std::unique_ptr<Bead<coordT>>& rhs) noexcept
{
    return distance_sq(lhs->position(), rhs->position());
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
min_distance(const std::unique_ptr<Bead<coordT>>& lhs,
             const std::unique_ptr<Bead<coordT>>& rhs) noexcept
{
    return std::sqrt(distance_sq(lhs, rhs));
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
min_distance_sq(const PDBAtom<coordT>& lhs, const PDBAtom<coordT>& rhs) noexcept
{
    return distance_sq(lhs->position(), rhs->position());
}

template<typename coordT>
inline typename scalar_type_of<coordT>::type
min_distance(const PDBAtom<coordT>& lhs, const PDBAtom<coordT>& rhs) noexcept
{
    return std::sqrt(distance_sq(lhs, rhs));
}

// containerT :=
//   {vector<PDBAtom>, vector<std::unique_ptr<Bead>>, PDBResidue, PDBChain}
// To add another class to this list, define ::value_type::coordinate_type.
template<typename containerT>
typename scalar_type_of<typename containerT::value_type::coordinate_type>::type
min_distance_sq(const containerT& lhs, const containerT& rhs)
{
    typedef typename scalar_type_of<
        typename containerT::value_type::coordinate_type>::type real_type;
    auto min_dist = std::numeric_limits<real_type>::infinity();

    for(auto l : lhs)
    {
        for(auto r : rhs)
        {
            const auto dist = min_distance_sq(l, r);
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

/* ---------------------------- distance_center ----------------------------- */

template<typename containerT>
typename scalar_type_of<typename containerT::value_type::coordinate_type>::type
distance_center_sq(const containerT& lhs, const containerT& rhs) noexcept
{
    return distance_sq(center(lhs), center(rhs));
}

template<typename containerT>
typename scalar_type_of<typename containerT::value_type::coordinate_type>::type
distance_center(const containerT& lhs, const containerT& rhs) noexcept
{
    return std::sqrt(distance_center_sq(lhs, rhs));
}

} // mjolnir
#endif /* JARNGREIPR_GEOMETRY_DISTANCE */
