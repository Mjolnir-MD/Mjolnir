#ifndef JARNGREIPR_IO_PDB_READER
#define JARNGREIPR_IO_PDB_READER
#include "PDBChain.hpp"
#include <fstream>
#include <sstream>

namespace jarngreipr
{

template<typename traitsT>
class PDBReader
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinatel_type coordinate_type;
    typedef PDBAtom<traits_type>    atom_type;
    typedef PDBResidue<traits_type> resiude_type;
    typedef PDBChain<traits_type>   chain_type;

   public:

    PDBReader()  = default;
    ~PDBReader() = default;

    std::vector<chain_type>
    parse(const std::vector<atom_type>& atoms) const;

    std::vector<atom_type>
    read(const std::string& fname) const;

    std::vector<atom_type>
    read(const std::string& fname, const std::size_t model_idx) const;

    std::vector<atom_type>
    read(std::basic_istream<char>& fname) const;

    std::vector<atom_type>
    read(std::basic_istream<char>& fname, const std::size_t model_idx) const;
};

template<typename traitsT>
std::vector<typename PDBReader<traitsT>::atom_type>
PDBReader<traitsT>::read(const std::string& fname) const
{
    std::ifstream ifs(fname);
    if(not ifs.good()) throw std::runtime_error("file open error: " + fname);
    const auto data = this->read(ifs);
    ifs.close();
    return data;
}

template<typename traitsT>
std::vector<typename PDBReader<traitsT>::atom_type>
PDBReader<traitsT>::read(
        const std::string& fname, const std::size_t model_idx) const
{
    std::ifstream ifs(fname);
    if(not ifs.good()) throw std::runtime_error("file open error: " + fname);
    const auto data = this->read(ifs, model_idx);
    ifs.close();
    return data;
}

template<typename traitsT>
std::vector<typename PDBReader<traitsT>::atom_type>
PDBReader<traitsT>::read(
        std::basic_istream<char>& ifs, const std::size_t model_idx) const
{
    while(!ifs.eof())
    {
        std::string line;
        std::getline(ifs, line);
        if(line.empty()) continue;
        if(line.substr(0, 5) == "MODEL")
        {
            std::size_t index;
            std::istringstream iss(line);
            std::string model;
            iss >> model >> index;
            if(index == model_idx) break;
        }
    }
    if(ifs.eof())
        throw std::invalid_argument("no model #" + std::to_string(model_idx));
    return this->read(ifs);
}

template<typename traitsT>
std::vector<typename PDBReader<traitsT>::atom_type>
PDBReader<traitsT>::read(std::basic_istream<char>& ifs) const
{
    std::vector<atom_type> atoms;
    while(not ifs.eof())
    {
        std::string line;
        std::getline(ifs, line);
        if(line.empty()) continue;
        if(line.substr(0, 3) == "END") break; // which is better ENDMDL or END?

        atom_type atom;
        if(line >> atom) atoms.emplace_back(std::move(atom));
    }
    return atoms;
}

template<typename traitsT>
std::vector<typename PDBReader<traitsT>::chain_type>
PDBReader<traitsT>::parse(const std::vector<atom_type>& atoms) const
{
    std::vector<residue_type> residues;
    residue_type tmp_residue;
    int current_residue_id = std::numeric_limits<int>::min();
    for(auto iter = atoms.cbegin(); iter != atoms.cend(); ++iter)
    {
        if(iter->residue_id != current_residue_id)
        {
            if(not tmp_residue.empty())
                residues.push_back(tmp_residue);
            tmp_residue.clear();
            current_residue_id = iter->residue_id;
        }
        tmp_residue.push_back(*iter);
    }
    if(not tmp_residue.empty()) residues.push_back(tmp_residue);

    std::vector<chain_type> chains;
    chain_type tmp_chain;
    std::string current_chain_id = "";
    for(auto iter = residues.begin(); iter != residues.end(); ++iter)
    {
        if(iter->chain_id() != current_chain_id)
        {
            if(not tmp_chain.empty()) chains.push_back(tmp_chain);
            tmp_chain.clear();
            current_chain_id = iter->chain_id();
        }
        tmp_chain.push_back(*iter);
    }
    if(not tmp_chain.empty()) chains.push_back(tmp_chain);

    return chains;
}

} // jarngreipr

#endif /* JARNGREIPR_IO_PDB_READER */
