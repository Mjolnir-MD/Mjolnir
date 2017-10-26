#ifndef JARNGREIPR_IO_PDB_WRITER
#define JARNGREIPR_IO_PDB_WRITER
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <jarngreipr/util/string.hpp>
#include <fstream>

namespace mjolnir
{

template<typename coordT>
std::ostream& operator<<(std::ostream& os, const PDBAtom<coordT>& atom)
{
    os << std::setw(6) << std::left << atom.prefix;
    os << std::setw(5) << std::right << atom.atom_id;
    os << ' ';
    if(atom.atom_name.size() < 4)
    {
        os << ' ' << std::setw(3) << std::left << atom.atom_name;
    }
    else // if(atom.atom_name.size() == 4)
    {
        os << std::setw(4) << atom.atom_name;
    }

    os << std::setw(1) << atom.altloc;
    os << std::setw(3) << std::right << atom.residue_name;
    os << ' ';
    os << std::setw(1) << atom.chain_id;
    os << std::setw(4) << std::right << atom.residue_id;
    os << std::setw(1) << atom.icode;
    os << "   ";
    os << std::setw(8) << std::fixed << std::setprecision(3) << std::right
       << atom.position[0];
    os << std::setw(8) << std::fixed << std::setprecision(3) << std::right
       << atom.position[1];
    os << std::setw(8) << std::fixed << std::setprecision(3) << std::right
       << atom.position[2];
    os << std::setw(6) << std::fixed << std::setprecision(2) << std::right
       << atom.occupancy;
    os << std::setw(6) << std::fixed << std::setprecision(2) << std::right
       << atom.temperature_factor;
    os << "          ";
    os << std::setw(2) << atom.element;
    os << std::setw(2) << atom.charge;
    return os;
}

template<typename coordT>
std::ostream& operator<<(std::ostream& os, const PDBResidue<coordT>& res)
{
    for(auto const& atom : res)
    {
        os << atom << '\n';
    }
    return os;
}

template<typename coordT>
std::ostream& operator<<(std::ostream& os, const PDBChain<coordT>& chain)
{
    for(auto const& res : chain)
    {
        os << res;
    }
    return os;
}

template<typename coordT>
void write_pdb_chains(const std::string& fname,
                      const std::vector<PDBChain<coordT>>& chains)
{
    std::ofstream ofs(fname);
    if(!ofs.good())
    {
        throw std::runtime_error(
                "mjolnir::io::read_pdb: file open error: " + fname);
    }
    for(const auto& chain : chains)
    {
        ofs << chain << '\n';
        ofs << "TER" << std::string(' ', 77);
    }
    ofs << std::flush;
    ofs.close();
    return;
}

template<typename coordT>
void write_pdb_atoms(const std::string& fname,
                     const std::vector<PDBAtom<coordT>>& atoms)
{
    std::ofstream ofs(fname);
    if(!ofs.good())
    {
        throw std::runtime_error(
                "mjolnir::io::read_pdb_atoms: file open error: " + fname);
    }
    for(const auto& atom : atoms)
    {
        ofs << atom << '\n';
    }
    ofs << std::flush;
    ofs.close();
    return;
}

} // mjolnir
#endif /* JARNGREIPR_IO_PDB_WRITER */
