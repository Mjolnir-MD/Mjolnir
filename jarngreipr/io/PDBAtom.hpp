#ifndef JARNGREIPR_IO_PDB_ATOM
#define JARNGREIPR_IO_PDB_ATOM
#include <jarngreipr/util/string.hpp>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdint>

namespace jarngreipr
{

template<typename traitsT>
struct PDBAtom
{
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::int64_t int_type;

    char            altloc;
    char            icode;
    int_type        atom_id;
    int_type        residue_id;
    real_type       occupancy;
    real_type       temperature_factor;
    std::string     prefix;
    std::string     atom_name;
    std::string     residue_name;
    std::string     chain_id;
    std::string     element;
    std::string     charge;
    coordinate_type position;

    explicit PDBAtom(const coordinate_type& pos);
    PDBAtom(const int_type aid, const coordinate_type& pos);
    PDBAtom(const int_type aid, const int_type resid, const std::string& atomname,
            const std::string& resname, const std::string& chainid,
            const coordinate_type& pos);

    PDBAtom() = default;
    ~PDBAtom() = default;
    PDBAtom(const PDBAtom&) = default;
    PDBAtom(PDBAtom&&)      = default;
    PDBAtom& operator=(const PDBAtom&) = default;
    PDBAtom& operator=(PDBAtom&&)      = default;
};

template<typename traitsT>
PDBAtom<traitsT>::PDBAtom(const coordinate_type& pos)
    : altloc(' '), icode(' '), atom_id(0), residue_id(0), occupancy(0.0),
      temperature_factor(0.0), prefix("ATOM"), atom_name("X"),
      residue_name("XXX"), chain_id("A"), element("X"), charge("X"),
      position(pos)
{
    // do nothing
}

template<typename traitsT>
PDBAtom<traitsT>::PDBAtom(const int_type aid, const coordinate_type& pos)
    : altloc(' '), icode(' '), atom_id(aid), residue_id(0), occupancy(0.0),
      temperature_factor(0.0), prefix("ATOM"), atom_name("X"),
      residue_name("XXX"), chain_id("A"), element("X"), charge("X"),
      position(pos)
{
    // do nothing
}

template<typename traitsT>
PDBAtom<traitsT>::PDBAtom(const int_type aid, const int_type resid, 
        const std::string& atomname, const std::string& resname,
        const std::string& chainid, const coordinate_type& pos)
    : altloc(' '), icode(' '), atom_id(aid), residue_id(resid), occupancy(0.0),
      temperature_factor(0.0), prefix("ATOM"), atom_name(atomname),
      residue_name(resname), chain_id(chainid), position(pos)
{
    // do nothing
}

template<typename traitsT>
std::basic_ostream<char>&
operator<<(std::basic_ostream<char>& os, const PDBAtom<traitsT>& a)
{
    os << std::setw(6) << std::left << a.prefix;
    os << std::setw(5) << std::right << a.atom_id;
    os << " ";
    os << std::setw(4) << a.atom_name;

    os << std::setw(1) << a.altloc;
    os << std::setw(3) << std::right << a.residue_name;
    os << " ";
    os << std::setw(1) << a.chain_id;
    os << std::setw(4) << std::right << a.residue_id;
    os << std::setw(1) << a.icode;
    os << "   ";
    os << std::setw(8) << std::fixed << std::setprecision(3) << std::right
       << a.position[0];
    os << std::setw(8) << std::fixed << std::setprecision(3) << std::right
       << a.position[1];
    os << std::setw(8) << std::fixed << std::setprecision(3) << std::right
       << a.position[2];
    os << std::setw(6) << std::fixed << std::setprecision(2) << std::right
       << a.occupancy;
    os << std::setw(6) << std::fixed << std::setprecision(2) << std::right
       << a.temperature_factor;
    os << "          ";
    os << std::setw(2) << a.element;
    os << std::setw(2) << a.charge;
    return os;
}

template<typename vector_type>
std::basic_istream<char>&
operator>>(std::basic_istream<char>& is, PDBAtom<vector_type>& atom)
{//XXX NOTE: read only one atom
    std::string line;
    while(!is.eof())
    {
        std::getline(is, line);
        if(line >> atom) return is;
        else continue;
    }
    return is;
}

template<typename traitsT>
bool operator>>(const std::string& line, PDBAtom<traitsT>& atom)
{
    typedef typename PDBAtom<traitsT>::real_type        real_type;
    typedef typename PDBAtom<traitsT>::coordinate_type coordinate_type;

    const std::string pref = remove_whitespace(line.substr(0, 6));
    if(pref != "ATOM" && pref != "HETATM") return false;
    atom.prefix       = pref;
    atom.atom_id      = std::stoi(line.substr(6, 5));
    atom.atom_name    = remove_whitespace(line.substr(12, 4));
    atom.altloc       = line[16];
    atom.residue_name = remove_whitespace(line.substr(17, 3));
    atom.chain_id     = line[21];
    atom.residue_id   = std::stoi(line.substr(22, 4));
    atom.icode        = line[26];
    real_type x, y, z;
    x = std::stod(line.substr(30, 8));
    y = std::stod(line.substr(38, 8));
    z = std::stod(line.substr(46, 8));
    atom.position = coordinate_type(x, y, z);

    try{atom.occupancy = std::stod(line.substr(54, 6));}
    catch(std::exception& excpt){atom.occupancy = 0e0;}

    try{atom.temperature_factor = std::stod(line.substr(60, 6));}
    catch(std::exception& excpt){atom.temperature_factor = 0e0;}

    try{atom.element = line.substr(76,2);}
    catch(std::exception& excpt){atom.element = "";}

    try{atom.charge = line.substr(78,2);}
    catch(std::exception& excpt){atom.charge = "";}

    return true;
}

}//jarngreipr
#endif // JARNGREIPR_PDB_ATOM
