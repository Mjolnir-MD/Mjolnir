#ifndef JARNGREIPR_GEOMETRY_CENTER
#define JARNGREIPR_GEOMETRY_CENTER
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <mjolnir/util/scalar_type_of.hpp>

namespace mjolnir
{

template<typename coordT>
inline coordT center(const std::unique_ptr<Bead<coordT>>& bead)
{
    return bead->position();
}

template<typename coordT>
inline coordT center(const PDBAtom<coordT>& atom) noexcept
{
    return atom.position;
}

template<typename containerT>
typename containerT::value_type::coordinate_type
center(const containerT& beads)
{
    typedef typename containerT::value_type::coordinate_type coord_type;
    typedef typename scalar_type_of<coord_type>::type real_type;
    coord_type pos(0,0,0);
    for(const auto& bead : beads)
    {
        pos += center(bead);
    }
    return pos / static_cast<real_type>(beads.size());
}

template<typename containerT, typename UnaryPredicate>
typename containerT::value_type::coordinate_type
center_if(const containerT& beads, UnaryPredicate&& satisfy)
{
    typedef typename containerT::value_type::coordinate_type coord_type;
    typedef typename scalar_type_of<coord_type>::type real_type;
    std::size_t N = 0;
    coord_type  pos(0,0,0);
    for(const auto& bead : beads)
    {
        if(satisfy(bead))
        {
            pos += bead.position();
            ++N;
        }
    }
    return pos / static_cast<real_type>(N);
}

template<typename Iterator>
inline typename std::iterator_traits<Iterator>::value_type::coordinate_type
center(Iterator first, const Iterator last) noexcept
{
    typedef typename std::iterator_traits<Iterator>::value_type::coordinate_type
            coordinate_type;
    typedef typename coordinate_type::real_type real_type;

    std::size_t num_p = 0;
    coordinate_type  pos(0,0,0);
    for(; first != last; ++first)
    {
        pos += center(*first);
        ++num_p;
    }
    return pos / static_cast<real_type>(num_p);
}

} // mjolnir
#endif// JARNGREIPR_GEOMETRY_CENTER
