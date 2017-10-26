#ifndef JARNGREIPR_GEOMETRY_CENTER
#define JARNGREIPR_GEOMETRY_CENTER
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <mjolnir/util/scalar_type_of.hpp>

namespace mjolnir
{

template<typename coordT>
inline coordT center(const PDBAtom<coordT>& atom)
{
    return atom.position;
}

template<typename coordT>
inline coordT center(const std::vector<PDBAtom<coordT>>& atoms)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    coordT pos(0,0,0);
    for(const auto& atom : atoms)
    {
        pos += atom.position;
    }
    return pos / static_cast<real_type>(atoms.size());
}

template<typename coordT>
inline coordT center(const PDBResidue<coordT>& residue)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    coordT pos(0,0,0);
    for(const auto& atom : residue)
    {
        pos += atom.position;
    }
    return pos / static_cast<real_type>(residue.size());
}

template<typename coordT>
inline coordT center(const PDBChain<coordT>& chain)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    std::size_t num_p = 0;
    coordT pos(0,0,0);
    for(const auto& residue : chain)
    {
        for(const auto& atom : residue)
        {
            pos += atom.position;
            ++num_p;
        }
    }
    return pos / static_cast<real_type>(num_p);
}

template<typename Iterator>
inline typename std::iterator_traits<Iterator>::value_type::coord_type
center(Iterator first, const Iterator last)
{
    typedef typename std::iterator_traits<Iterator>::value_type::coord_type
            coord_type;
    typedef typename coord_type::real_type real_type;

    std::size_t N = 0;
    coord_type  pos(0,0,0);
    for(; first != last; ++first)
    {
        pos += center(*first);
        ++N;
    }
    return pos / static_cast<real_type>(N);
}

} // mjolnir
#endif// JARNGREIPR_GEOMETRY_CENTER
