#ifndef JARNGREIPR_IO_PDB_WRITER
#define JARNGREIPR_IO_PDB_WRITER
#include "PDBChain.hpp"
#include <fstream>

namespace jarngreipr
{

template<typename traitsT>
class PDBWriter
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinatel_type coordinate_type;
    typedef PDBAtom<traits_type>    atom_type;
    typedef PDBResidue<traits_type> resiude_type;
    typedef PDBChain<traits_type>   chain_type;

  public:

    PDBWriter()  = default;
    ~PDBWriter() = default;

    void write(const std::string& filename,
               const std::vector<atom_type>& atoms) const;

    void write(std::ostream& os, const std::vector<atom_type>& atoms) const;

    void write(const std::string& fname,
               const std::vector<chain_type>& model) const;

    void write(std::ostream& os, const std::vector<chain_type>& model) const;
};

template<typename traitsT>
void PDBWriter<traitsT>::write(
        const std::string& filename, const std::vector<atom_type>& atoms) const
{
    std::ofstream ofs(filename);
    if(not ofs.good()) throw std::runtime_error("file open error: " + filename);
    this->write(ofs, atoms);
    ofs.close();
    return;
}

template<typename traitsT>
void PDBWriter<traitsT>::write(
        const std::string& filename, const std::vector<chain_type>& model) const
{
    std::ofstream ofs(filename);
    if(not ofs.good()) throw std::runtime_error("file open error: " + filename);
    this->write(ofs, model);
    ofs.close();
    return;
}

template<typename traitsT>
void PDBWriter<traitsT>::write(
        std::ostream& os, const std::vector<atom_type>& atoms) const
{
    for(auto iter = atoms.cbegin(); iter != atoms.cend(); ++iter)
        os << *iter << std::endl;
    return;
}

template<typename traitsT>
void PDBWriter<traitsT>::write(
        std::ostream& os, const std::vector<chain_type>& chains) const
{
    std::size_t index = 1;
    for(auto chain = chains.cbegin(); chain != chains.cend(); ++chain)
    {
        for(auto resi = chain->cbegin(); resi != chain->cend(); ++resi)
            for(auto iter = resi->cbegin(); iter != resi->cend(); ++iter)
                os << *iter << std::endl;
        os << "TER" << std::endl;
        ++index;
    }
    return;
}

} // jarn
#endif /* JARNGREIPR_IO_PDB_WRITER */
