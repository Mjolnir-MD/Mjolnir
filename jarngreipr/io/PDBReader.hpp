#ifndef JARNGREIPR_IO_PDB_READER
#define JARNGREIPR_IO_PDB_READER
#include <jarngreipr/io/PDBAtom.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/io/PDBChain.hpp>
#include <jarngreipr/util/string.hpp>
#include <stdexcept>
#include <fstream>
#include <cassert>

namespace mjolnir
{

namespace pdb
{
struct different_line : public std::invalid_argument {
    explicit different_line(const std::string& w) : std::invalid_argument(w) {}
    explicit different_line(const char* w)        : std::invalid_argument(w) {}
    ~different_line() override = default;
};

struct invalid_format : public std::runtime_error {
    explicit invalid_format(const std::string& w) : std::runtime_error(w) {}
    explicit invalid_format(const char* w)        : std::runtime_error(w) {}
    ~invalid_format() override = default;
};

struct no_more_atom : public std::runtime_error {
    explicit no_more_atom(const std::string& w) : std::runtime_error(w) {}
    explicit no_more_atom(const char* w)        : std::runtime_error(w) {}
    ~no_more_atom() override = default;
};
}// pdb

template<typename coordT>
PDBAtom<coordT> read_pdb_atom_line(const std::string& line)
{
    try
    {
        const std::string prefix = remove_whitespace(line.substr(0, 6));
        if(prefix != "ATOM")
        {
            throw pdb::different_line(
                    "mjolnir::read_pdb_atom: neither ATOM nor HETATM line");
        }
        PDBAtom<coordT> atom;
        atom.prefix       = prefix;
        atom.atom_id      = std::stoi(line.substr(6, 5));
        atom.atom_name    = remove_whitespace(line.substr(12, 4));
        atom.altloc       = line[16];
        atom.residue_name = remove_whitespace(line.substr(17, 3));
        atom.chain_id     = line[21];
        atom.residue_id   = std::stoi(line.substr(22, 4));
        atom.icode        = line[26];
        typename PDBAtom<coordT>::real_type x, y, z;
        x = std::stod(line.substr(30, 8));
        y = std::stod(line.substr(38, 8));
        z = std::stod(line.substr(46, 8));
        atom.position = coordT(x, y, z);

        try{atom.occupancy = std::stod(line.substr(54, 6));}
        catch(std::exception& excpt){atom.occupancy = 0e0;}

        try{atom.temperature_factor = std::stod(line.substr(60, 6));}
        catch(std::exception& excpt){atom.temperature_factor = 0e0;}

        try{atom.element = line.substr(76,2);}
        catch(std::exception& excpt){atom.element = "";}

        try{atom.charge = line.substr(78,2);}
        catch(std::exception& excpt){atom.charge = "";}

        return atom;
    }
    catch(pdb::different_line const&)
    {
        throw;
    }
    catch(...)
    {
        throw pdb::invalid_format(
                "mjolnir::read_pdb_atom_line: reading line"_str + line);
    }
}

template<typename coordT>
std::vector<PDBAtom<coordT>> read_pdb_atoms(std::istream& istrm)
{
    std::vector<PDBAtom<coordT>> atoms;
    while(!istrm.eof())
    {
        std::string line;
        std::getline(istrm, line);
        try
        {
            atoms.push_back(read_pdb_atom_line<coordT>(line));
        }
        catch(pdb::different_line const&)
        {
            // found non-ATOM line. ignore and continue.
            continue;
        }
        catch(pdb::invalid_format const& ivf)
        {
            throw std::runtime_error("mjolnir::io::read_pdb_atoms: "_str +
                "invalid line: \n"_str + line);
        }
    }
    return atoms;
}

template<typename coordT>
std::vector<PDBChain<coordT>> read_pdb_chains(std::istream& istrm)
{
    const std::vector<PDBAtom<coordT>> atoms = read_pdb_atoms<coordT>(istrm);

    std::vector<PDBResidue<coordT>> residues;
    /* split-atoms-into-residues */{
        PDBResidue<coordT> tmp;
        for(const auto& atom : atoms)
        {
            if(tmp.is_in_same_residue(atom))
            {
                tmp.push_back(atom);
            }
            else // make new residue
            {
                assert(false == tmp.empty());

                residues.push_back(tmp);
                tmp.clear();
                tmp.push_back(atom);
            }
        }
        if(false == tmp.empty())
        {
            residues.push_back(std::move(tmp));
        }
    }

    std::vector<PDBChain<coordT>> chains;
    /* split-residues-into-chians */{
        PDBChain<coordT> tmp;
        for(auto&& residue : residues)
        {
            if(tmp.is_in_same_chain(residue))
            {
                tmp.push_back(std::move(residue));
            }
            else // make new chain
            {
                assert(false == tmp.empty());

                chains.push_back(tmp);
                tmp.clear();
                tmp.push_back(residue);
            }
        }
        if(false == tmp.empty())
        {
            chains.push_back(std::move(tmp));
        }
    }
    return chains;
}

template<typename coordT>
std::vector<PDBChain<coordT>> read_pdb_chains(const std::string& fname)
{
    std::ifstream ifs(fname);
    if(!ifs.good())
    {
        throw std::runtime_error(
                "mjolnir::io::read_pdb_chains: file open error: " + fname);
    }
    return read_pdb_chains<coordT>(ifs);
}

template<typename coordT>
std::vector<PDBAtom<coordT>> read_pdb_atoms(const std::string& fname)
{
    std::ifstream ifs(fname);
    if(!ifs.good())
    {
        throw std::runtime_error(
                "mjolnir::io::read_pdb_atoms: file open error: " + fname);
    }
    return read_pdb_atoms<coordT>(ifs);
}


// lazy PDB reader
template<typename coordT>
class PDBReader
{
  public:

    typedef coordT coordinate_type;
    typedef PDBAtom<coordinate_type>    atom_type;
    typedef PDBResidue<coordinate_type> residue_type;
    typedef PDBChain<coordinate_type>   chain_type;
  public:

    explicit PDBReader(const std::string& fname)
        : filename_(fname), ifs_(fname)
    {}
    ~PDBReader() = default;

    atom_type    read_next_atom();
    residue_type read_next_residue();
    chain_type   read_next_chain();

    std::string const& filename() const noexcept {return this->filename_;}

    bool is_eof()  {ifs_.peek(); return ifs_.eof();}
    bool is_good() {return ifs_.good();}

  private:

    const std::string filename_;
    std::ifstream     ifs_;
};

template<typename coordT>
PDBAtom<coordT> PDBReader<coordT>::read_next_atom()
{
    while(!this->ifs_.eof())
    {
        std::string line;
        std::getline(this->ifs_, line);
        if(line.size() < 6)
        {
            // cannot detect ATOM or HERATM signeture. ignore the line.
            continue;
        }
        try
        {
            return read_pdb_atom_line<coordT>(line);
        }
        catch(pdb::different_line const&)
        {
            continue;
        }
        catch(pdb::invalid_format const&)
        {
            throw pdb::invalid_format(
                    "mjolnir::io::PDBReader::read_next_atom: "_str +
                    "invalid line found in the file: \n"_str + line);
        }
    }
    throw pdb::no_more_atom("mjolnir::PDBReader::read_next_atom: "_str +
            "no more atoms in file: "_str + this->filename_);
}

template<typename coordT>
PDBResidue<coordT> PDBReader<coordT>::read_next_residue()
{
    residue_type res;
    while(!this->ifs_.eof())
    {
        const auto pos = this->ifs_.tellg();
        atom_type atm;
        try
        {
            atm = this->read_next_atom();
        }
        catch(const pdb::no_more_atom& nma)
        {
            break;
        }

        if(res.is_in_same_residue(atm))
        {
            res.push_back(atm);
        }
        else // different residue
        {
            this->ifs_.seekg(pos); // restore the position
            break;
        }
    }
    if(res.empty())
    {
        throw pdb::no_more_atom("mjolnir::PDBReader::read_next_residue: "_str +
                "no more atoms in file: "_str + this->filename_);
    }
    return res;
}

template<typename coordT>
PDBChain<coordT> PDBReader<coordT>::read_next_chain()
{
    chain_type chn;
    while(!this->ifs_.eof())
    {
        const auto pos = this->ifs_.tellg();
        residue_type res;
        try
        {
            res = this->read_next_residue();
        }
        catch(const pdb::no_more_atom& nma)
        {
            break;
        }

        if(chn.is_in_same_chain(res))
        {
            chn.push_back(res);
        }
        else // different residue
        {
            this->ifs_.seekg(pos); // restore the position
            break;
        }
    }
    if(chn.empty())
    {
        throw pdb::no_more_atom("mjolnir::PDBReader::read_next_chain: "_str +
                "no more atoms in file: "_str + this->filename_);
    }
    return chn;
}

} // mjolnir
#endif /* JARNGREIPR_IO_PDB_READER */
