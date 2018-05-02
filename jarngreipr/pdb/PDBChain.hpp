#ifndef JARNGREIPR_IO_PDB_CHAIN
#define JARNGREIPR_IO_PDB_CHAIN
#include <jarngreipr/pdb/PDBAtom.hpp>
#include <mjolnir/util/range.hpp>
#include <vector>

namespace jarngreipr
{

template<typename realT, typename coordT>
class PDBChain
{
  public:
    typedef PDBAtom<realT, coordT> atom_type;
    typedef std::vector<atom_type> container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

    typedef mjolnir::range<const_iterator>   const_residue_range;
    typedef std::vector<const_residue_range> residue_container;
    typedef typename residue_container::iterator       residue_iterator;
    typedef typename residue_container::const_iterator const_residue_iterator;

  public:

    explicit PDBChain(std::vector<atom_type> atoms): atoms_(std::move(atoms))
    {
        std::int32_t residue_id = this->atoms_.front().residue_id;
        const_iterator    first = this->atoms_.begin();
        for(auto i = this->atoms_.begin(), e = this->atoms_.end(); i!=e; ++i)
        {
            if(i->residue_id != residue_id)
            {
                this->residues_.push_back(const_residue_range(first, i));
                first = i;
                residue_id = i->residue_id;
            }
        }
        this->residues_.push_back(const_residue_range(first, this->atoms_.end()));
    }

    ~PDBChain() = default;
    PDBChain(const PDBChain&) = default;
    PDBChain(PDBChain&&)      = default;
    PDBChain& operator=(const PDBChain&) = default;
    PDBChain& operator=(PDBChain&&)      = default;

    char chain_id() const noexcept {return this->atoms_.front().chain_id;}

    bool empty() const noexcept {return atoms_.empty();}
    void clear()                {return atoms_.clear();}
    void resize (const std::size_t s){atoms_.resize(s);}
    void reserve(const std::size_t s){atoms_.reserve(s);}

    std::size_t atoms_size()    const noexcept {return atoms_.size();}
    std::size_t residues_size() const noexcept {return residues_.size();}

    const_residue_range residue_at(const std::size_t i) const {return residues_.at(i);}
    const_residue_range atom_at   (const std::size_t i) const {return atoms_.at(i);}

    const_iterator begin()  const noexcept {return atoms_.begin();}
    const_iterator end()    const noexcept {return atoms_.end();}
    const_iterator cbegin() const noexcept {return atoms_.cbegin();}
    const_iterator cend()   const noexcept {return atoms_.cend();}

    const_residue_iterator res_begin()  const noexcept {return residues_.begin();}
    const_residue_iterator res_end()    const noexcept {return residues_.end();}
    const_residue_iterator res_cbegin() const noexcept {return residues_.cbegin();}
    const_residue_iterator res_cend()   const noexcept {return residues_.cend();}

  private:

    std::vector<const_residue_range> residues_;
    container_type atoms_;
};

}//mjolnir
#endif /* JARNGREIPR_IO_PDB_CHAIN */
