#ifndef JARNGREIPR_PDB_READER_HPP
#define JARNGREIPR_PDB_READER_HPP
#include <jarngreipr/pdb/PDBAtom.hpp>
#include <jarngreipr/pdb/PDBChain.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <fstream>
#include <sstream>

namespace jarngreipr
{

// lazy pdb reader
template<typename realT>
class PDBReader
{
  public:
    typedef PDBAtom<realT>  atom_type;
    typedef PDBChain<realT> chain_type;
    typedef std::vector<chain_type> model_type;

  public:

    explicit PDBReader(const std::string& fname)
        : filename_(fname), ifstrm_(fname)
    {
        if(!ifstrm_.good())
        {
            throw std::runtime_error(
                "jarngreipr::PDBReader: file open error: " + fname);
        }
    }

    bool is_eof() const noexcept {return this->ifstrm_.eof();}

    void rewind() {this->ifstrm_.seekg(0, std::ios::beg);}

    chain_type read_chain(const char id)
    {
        this->rewind();
        while(!this->ifstrm_.eof())
        {
            try
            {
                const auto chain = this->read_next_chain();
                if(chain.chain_id() == id)
                {
                    return chain;
                }
            }
            catch(std::runtime_error const& re)
            {
                break;
            }
        }
        mjolnir::throw_exception<std::runtime_error>("file ",
                this->filename_, " does not contain chain ", id);
    }

    // lazy functions. throws std::runtime_error if it reaches EOF.
    atom_type read_next_atom()
    {
        atom_type atm;
        while(!this->ifstrm_.eof())
        {
            std::string line;
            std::getline(ifstrm_, line);
            ifstrm_.peek(); // set eof flag if it reached

            std::istringstream iss(line);
            try
            {
                iss >> atm;
            }
            catch(std::runtime_error const& re)
            {
                continue;
            }
            return atm;
        }
        mjolnir::throw_exception<std::runtime_error>("file ", this->filename_,
            " does not contain ATOM any more.");
    }
    chain_type read_next_chain()
    {
        std::vector<atom_type> atoms;
        while(!this->ifstrm_.eof())
        {
            atom_type atm;
            std::string line;
            std::getline(ifstrm_, line);
            ifstrm_.peek(); // set eof flag

            std::istringstream iss(line);
            try
            {
                iss >> atm;
                atoms.push_back(atm);
            }
            catch(std::runtime_error const& re)
            {
                if(line.substr(0, 3) == "TER" || line.substr(0, 6) == "ENDMDL")
                {
                    return chain_type(std::move(atoms));
                }
            }
            continue;
        }
        if(!atoms.empty())
        {
            return chain_type(std::move(atoms));
        }
        mjolnir::throw_exception<std::runtime_error>("file ", this->filename_,
            " does not contain chain any more.");
    }

    model_type read_next_model()
    {
        std::vector<chain_type> chains;
        std::vector<atom_type>  atoms;
        while(!this->ifstrm_.eof())
        {
            atom_type atm;
            std::string line;
            std::getline(ifstrm_, line);
            ifstrm_.peek(); // set eof flag

            std::istringstream iss(line);
            try
            {
                iss >> atm;
                atoms.push_back(atm);
            }
            catch(std::runtime_error const& re)
            {
                if(line.substr(0, 3) == "TER" || line.substr(0, 6) == "ENDMDL")
                {
                    chains.push_back(chain_type(std::move(atoms)));
                    atoms = {};
                }
                if(line.substr(0, 6) == "ENDMDL")
                {
                    return chains;
                }
            }
            continue;
        }
        if(!chains.empty())
        {
            return chains;
        }
        mjolnir::throw_exception<std::runtime_error>("file ", this->filename_,
            " does not contain chain any more.");
    }

  private:
    std::string filename_;
    std::ifstream ifstrm_;
};

} // jarngreipr
#endif// JARNGREIPR_PDB_READER_HPP
