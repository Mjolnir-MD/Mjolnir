#ifndef JARNGREIPR_PDB_ATOM_HPP
#define JARNGREIPR_PDB_ATOM_HPP
#include <mjolnir/math/Vector.hpp>
#include <ostream>
#include <istream>
#include <iomanip>
#include <string>
#include <cstdint>

namespace jarngreipr
{

template<typename realT>
struct PDBAtom
{
    using real_type   = realT;
    using coordinate_type = mjolnir::Vector<realT, 3>;

    char         altloc;
    char         icode;
    char         chain_id;
    std::int32_t atom_id;
    std::int32_t residue_id;
    real_type    occupancy;
    real_type    temperature_factor;
    std::string  atom_name; //XXX with whitespaces!
    std::string  residue_name;
    std::string  element;
    std::string  charge;
    coordinate_type position;
};

template<typename realT>
PDBAtom<realT> make_pdb_atom(const mjolnir::Vector<realT, 3>& pos,
    const std::int32_t atm_id =   0, const std::string& atm_name = " X  ",
    const std::int32_t res_id =   0, const std::string& res_name = "XXX",
    const char         chn_id = 'A', const std::string& elm_name = "X")
{
    return PDBAtom<realT>{' ', ' ', chn_id, atm_id, res_id, 0.0, 0.0,
                          atm_name, res_name, elm_name, "", pos};
}

template<typename charT, typename traits, typename realT>
std::basic_ostream<charT, traits>& operator<<(
    std::basic_ostream<charT, traits>& os, const PDBAtom<realT>& atm)
{
    os << "ATOM  ";
    os << std::right << std::setw(5) << atm.atom_id;
    os << ' ';
    os << std::setw(4) << atm.atom_name;
    os << atm.altloc;
    os << std::left  << std::setw(3) << atm.residue_name;
    os << ' ';
    os << atm.chain_id;
    os << std::right << std::setw(4) << atm.residue_id;
    os << atm.icode;
    os << "   ";
    os << std::right << std::setw(8) << std::fixed << std::setprecision(3)
       << atm.position[0];
    os << std::right << std::setw(8) << std::fixed << std::setprecision(3)
       << atm.position[1];
    os << std::right << std::setw(8) << std::fixed << std::setprecision(3)
       << atm.position[2];
    os << std::right << std::setw(6) << std::fixed << std::setprecision(2)
       << atm.occupancy;
    os << std::right << std::setw(6) << std::fixed << std::setprecision(2)
       << atm.temperature_factor;
    os << "          ";
    os << std::right << std::setw(2) << atm.element;
    os << std::right << std::setw(2) << atm.charge;
    return os;
}

template<typename charT, typename traits, typename realT>
std::basic_istream<charT, traits>& operator>>(
    std::basic_istream<charT, traits>& is, PDBAtom<realT>& atm)
{
    std::string line;
    std::getline(is, line);
    if(line.substr(0, 6) != "ATOM  ")
    {
        throw std::runtime_error("not an ATOM line: " + line);
    }
    atm.atom_id      = std::stoi(line.substr( 6, 5));
    atm.atom_name    = line.substr(12, 4);
    atm.altloc       = line.at(16);
    atm.residue_name = line.substr(17, 3);
    atm.chain_id     = line.at(21);
    atm.residue_id   = std::stoi(line.substr(22, 4));
    atm.icode        = line.at(26);
    atm.position[0]  = std::stod(line.substr(30, 8));
    atm.position[1]  = std::stod(line.substr(38, 8));
    atm.position[2]  = std::stod(line.substr(46, 8));

    atm.occupancy          = 0.0;
    atm.temperature_factor = 0.0;
    atm.element            = "  ";
    atm.charge             = "  ";

    try{atm.occupancy = std::stod(line.substr(54, 6));}
    catch(const std::out_of_range& o){return is;}
    try{atm.temperature_factor = std::stod(line.substr(60, 6));}
    catch(const std::out_of_range& o){return is;}
    try{atm.element = line.substr(76, 2);}
    catch(const std::out_of_range& o){return is;}
    try{atm.charge = line.substr(78, 2);}
    catch(const std::out_of_range& o){return is;}
    return is;
}

}//jarngreipr
#endif // JARNGREIPR_PDB_ATOM_HPP
