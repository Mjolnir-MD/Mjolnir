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
inline coordT center(const Bead<coordT>& bead)
{
    return bead.position();
}

template<typename coordT>
coordT center(const std::vector<Bead<coordT>>& beads)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    coordT pos(0,0,0);
    for(const auto& bead : beads)
    {
        pos += bead.position();
    }
    return pos / static_cast<real_type>(beads.size());
}

template<typename coordT, typename UnaryPredicate>
coordT center_if(const std::vector<Bead<coordT>>& beads,
                 UnaryPredicate&& satisfy)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    std::size_t N = 0;
    coordT pos(0,0,0);
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


template<typename coordT>
inline coordT center(const PDBAtom<coordT>& atom) noexcept
{
    return atom.position;
}

template<typename coordT>
coordT center(const std::vector<PDBAtom<coordT>>& atoms) noexcept
{
    typedef typename scalar_type_of<coordT>::type real_type;
    coordT pos(0,0,0);
    for(const auto& atom : atoms)
    {
        pos += atom.position;
    }
    return pos / static_cast<real_type>(atoms.size());
}

template<typename coordT, typename UnaryPredicate>
coordT center_if(const std::vector<PDBAtom<coordT>>& atoms,
                 UnaryPredicate&& satisfy)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    std::size_t N = 0;
    coordT pos(0,0,0);
    for(const auto& atom : atoms)
    {
        if(satisfy(atom))
        {
            pos += atom.position;
            ++N;
        }
    }
    return pos / static_cast<real_type>(N);
}

template<typename coordT>
coordT center(const PDBResidue<coordT>& residue) noexcept
{
    typedef typename scalar_type_of<coordT>::type real_type;
    coordT pos(0,0,0);
    for(const auto& atom : residue)
    {
        pos += atom.position;
    }
    return pos / static_cast<real_type>(residue.size());
}

template<typename coordT, typename UnaryPredicate>
coordT center_if(const PDBResidue<coordT>& residue, UnaryPredicate&& satisfy)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    std::size_t N = 0;
    coordT pos(0,0,0);
    for(const auto& atom : residue)
    {
        if(satisfy(atom))
        {
            pos += atom.position;
            ++N;
        }
    }
    return pos / static_cast<real_type>(N);
}


template<typename coordT>
coordT center(const PDBChain<coordT>& chain) noexcept
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

template<typename coordT, typename UnaryPredicate>
coordT center_if(const PDBChain<coordT>& chain, UnaryPredicate&& satisfy)
{
    typedef typename scalar_type_of<coordT>::type real_type;
    std::size_t num_p = 0;
    coordT pos(0,0,0);
    for(const auto& residue : chain)
    {
        for(const auto& atom : residue)
        {
            if(satisfy(atom))
            {
                pos += atom.position;
                ++num_p;
            }
        }
    }
    return pos / static_cast<real_type>(num_p);
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
