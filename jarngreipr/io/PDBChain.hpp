#ifndef JARNGREIPR_IO_PDB_CHAIN
#define JARNGREIPR_IO_PDB_CHAIN
#include <jarngreipr/io/PDBResidue.hpp>

namespace mjolnir
{

template<typename coordT>
class PDBChain
{
  public:
    typedef PDBResidue<coordT> residue_type;
    typedef typename residue_type::int_type        int_type;
    typedef typename residue_type::real_type       real_type;
    typedef typename residue_type::coordinate_type coordinate_type;
    typedef typename residue_type::atom_type       atom_type;

    typedef std::vector<residue_type> container_type;
    typedef typename container_type::value_type     value_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    PDBChain()  = default;
    ~PDBChain() = default;
    PDBChain(const PDBChain&) = default;
    PDBChain(PDBChain&&)      = default;
    PDBChain& operator=(const PDBChain&) = default;
    PDBChain& operator=(PDBChain&&)      = default;

    char chain_id() const noexcept {return this->front().chain_id();}

    bool        empty() const noexcept {return residues_.empty();}
    std::size_t size()  const noexcept {return residues_.size();}
    void clear(){return residues_.clear();}
    void resize (const std::size_t s){residues_.resize(s);}
    void reserve(const std::size_t s){residues_.reserve(s);}

    void push_back(const residue_type& res)
    {
        if(false == is_in_same_chain(res))
        {
            throw std::invalid_argument("mjolnir::PDBChain::push_back:"
                    "pushing a residue in different chain");
        }
        residues_.push_back(res);
        return;
    }
    void push_back(residue_type&& res)
    {
        if(false == is_in_same_chain(res))
        {
            throw std::invalid_argument("mjolnir::PDBChain::push_back:"
                    "pushing a residue in different chain");
        }
        residues_.push_back(std::move(res));
        return;
    }

    residue_type const& at(std::size_t i) const {return residues_.at(i);}
    residue_type &      at(std::size_t i)       {return residues_.at(i);}
    residue_type const& operator[](std::size_t i) const noexcept
    {return residues_[i];}
    residue_type &      operator[](std::size_t i)       noexcept
    {return residues_[i];}

    residue_type const& front() const noexcept {return residues_.front();}
    residue_type &      front()       noexcept {return residues_.front();}
    residue_type const& back()  const noexcept {return residues_.back();}
    residue_type &      back()        noexcept {return residues_.back();}

    iterator       begin()        noexcept {return residues_.begin();}
    iterator       end()          noexcept {return residues_.end();}
    const_iterator begin()  const noexcept {return residues_.begin();}
    const_iterator end()    const noexcept {return residues_.end();}
    const_iterator cbegin() const noexcept {return residues_.cbegin();}
    const_iterator cend()   const noexcept {return residues_.cend();}

    container_type&       residues()       noexcept {return residues_;}
    container_type const& residues() const noexcept {return residues_;}

    bool is_in_same_chain(const residue_type& r) const noexcept
    {
        return this->empty() || (r.chain_id() == this->chain_id());
    }

  private:

    container_type residues_;
};

}//mjolnir
#endif /* JARNGREIPR_IO_PDB_CHAIN */
