#ifndef JARNGREIPR_PDB_WRITER_HPP
#define JARNGREIPR_PDB_WRITER_HPP
#include <jarngreipr/pdb/PDBAtom.hpp>
#include <jarngreipr/pdb/PDBChain.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <fstream>
#include <sstream>

namespace jarngreipr
{

template<typename realT>
class PDBWriter
{
  public:
    typedef PDBAtom<realT>  atom_type;
    typedef PDBChain<realT> chain_type;
    typedef std::vector<chain_type> model_type;

  public:

    explicit PDBWriter(const std::string& fname): filename_(fname), ofstrm_(fname)
    {
        if(!ofstrm_.good())
        {
            throw std::runtime_error(
                "jarngreipr::PDBWriter: file open error: " + fname);
        }
    }

    void write_chain(const chain_type& chain)
    {
        for(const auto& atom : chain)
        {
            ofstrm_ << atom << '\n';
        }
        ofstrm_ << "TER\n";
    }

    void write_model(const std::vector<chain_type>& model,
                     const std::size_t N_model)
    {
        ofstrm_ << "MODEL " << N_model << '\n';
        for(const auto& chain : model)
        {
            this->write_chain(chain);
        }
        ofstrm_ << "ENDMDL\n";
    }

  private:
    std::string filename_;
    std::ofstream ofstrm_;
};

} // jarngreipr
#endif// JARNGREIPR_PDB_WRITER_HPP
