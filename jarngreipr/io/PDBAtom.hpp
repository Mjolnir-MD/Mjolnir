#ifndef JARNGREIPR_IO_PDB_ATOM
#define JARNGREIPR_IO_PDB_ATOM
#include <mjolnir/util/scalar_type_of.hpp>
#include <jarngreipr/util/string.hpp>
#include <type_traits>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdint>

namespace mjolnir
{

template<typename coordT>
struct PDBAtom
{
    typedef coordT   coord_type;
    typedef typename scalar_type_of<coord_type>::type scalar_type;
    typedef scalar_type real_type;
    typedef std::int64_t int_type;

    static_assert(std::is_floating_point<scalar_type>::value,
                  "scalar_type_of<PDBAtom::coordT> should be floating point");

    char        altloc;
    char        icode;
    int_type    atom_id;
    int_type    residue_id;
    real_type   occupancy;
    real_type   temperature_factor;
    std::string prefix;
    std::string atom_name;
    std::string residue_name;
    std::string chain_id;
    std::string element;
    std::string charge;
    coord_type  position;

    explicit PDBAtom(const coord_type& pos);
    PDBAtom(const int_type atomid, const coord_type& pos);
    PDBAtom(const int_type atomid, const int_type resid,
            const std::string& atomname, const std::string& resname,
            const std::string& chainid,  const coord_type& pos);

    PDBAtom()  = default;
    ~PDBAtom() = default;
    PDBAtom(const PDBAtom&) = default;
    PDBAtom(PDBAtom&&)      = default;
    PDBAtom& operator=(const PDBAtom&) = default;
    PDBAtom& operator=(PDBAtom&&)      = default;
};

template<typename coordT>
PDBAtom<coordT>::PDBAtom(const coord_type& pos)
    : altloc(' '), icode(' '), atom_id(0), residue_id(0), occupancy(0.0),
      temperature_factor(0.0), prefix("ATOM"), atom_name("X"),
      residue_name("XXX"), chain_id("X"), element("X"), charge("X"),
      position(pos)
{}

template<typename coordT>
PDBAtom<coordT>::PDBAtom(const int_type atomid, const coord_type& pos)
    : altloc(' '), icode(' '), atom_id(atomid), residue_id(0), occupancy(0.0),
      temperature_factor(0.0), prefix("ATOM"), atom_name("X"),
      residue_name("XXX"), chain_id("A"), element("X"), charge("X"),
      position(pos)
{}

template<typename coordT>
PDBAtom<coordT>::PDBAtom(const int_type atomid, const int_type resid,
        const std::string& atomname, const std::string& resname,
        const std::string& chainid,  const coord_type& pos)
    : altloc(' '), icode(' '), atom_id(atomid), residue_id(resid),
      occupancy(0.0), temperature_factor(0.0), prefix("ATOM"),
      atom_name(atomname), residue_name(resname), chain_id(chainid),
      position(pos)
{}

}//jarngreipr
#endif // JARNGREIPR_PDB_ATOM
