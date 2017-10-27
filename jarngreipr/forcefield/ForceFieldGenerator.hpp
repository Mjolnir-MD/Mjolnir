#ifndef JARNGREIPR_FORCEFIELD_GENERATOR
#define JARNGREIPR_FORCEFIELD_GENERATOR
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/model/CGChain.hpp>
#include <jarngreipr/model/CarbonAlpha.hpp>
#include <jarngreipr/forcefield/ConnectingIndices.hpp>
#include <ostream>

namespace mjolnir
{

template<typename coordT>
class ForceFieldGenerator
{
  public:
    typedef coordT                      coordinate_type;
    typedef PDBAtom<coordinate_type>    atom_type;
    typedef PDBResidue<coordinate_type> residue_type;
    typedef PDBChain<coordinate_type>   chain_type;
    typedef Bead<coordinate_type>       bead_type;
    typedef CGChain<coordinate_type>    cg_chain_type;
    typedef std::map<std::size_t, ConnectingIndices> connection_info;

  public:
    virtual ~ForceFieldGenerator() = default;

    //!@brief generate parameter values and write out to ostream.
    //@return return connection info. for scaling or ignoring Lennard-Jones.
    virtual connection_info
    generate(std::ostream& ostrm, const cg_chain_type& chain) const = 0;

    //!@brief if chain contains invalid bead, return false.
    virtual bool
    check_bead_kind(const cg_chain_type& chain) const = 0;
};

} // mjolnir
#endif// JARNGREIPR_FORCEFIELD_GENERATOR
