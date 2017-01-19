#ifndef JARNGREIPR_IO_PDB_CHAIN
#define JARNGREIPR_IO_PDB_CHAIN
#include "PDBResidue.hpp"

namespace jarngreipr
{

template<typename traitsT>
class PDBChain
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef PDBResidue<traits_type> residue_type;
    typedef typename residue_type::int_type int_type;
    typedef std::vector<residue_type> container_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    PDBChain()  = default;
    ~PDBChain() = default;
    PDBChain(const PDBChain&) = default;
    PDBChain(PDBChain&&)      = default;
    PDBChain& operator=(const PDBChain&) = default;
    PDBChain& operator=(PDBChain&&)      = default;

    std::string const& chain_id() const {return this->front().front().chain_id;}

    bool empty() const {return residues_.empty();}
    std::size_t size() const {return residues_.size();}
    void push_back(const residue_type& res);
    void emplace_back(residue_type&& res);
    void clear(){return residues_.clear();}

    residue_type const& at(std::size_t i) const {return residues_.at(i);}
    residue_type &      at(std::size_t i)       {return residues_.at(i);}
    residue_type const& operator[](std::size_t i) const {return residues_[i];}
    residue_type &      operator[](std::size_t i)       {return residues_[i];}

    residue_type const& front() const {return residues_.front();}
    residue_type &      front()       {return residues_.front();}
    residue_type const& back() const {return residues_.back();}
    residue_type &      back()       {return residues_.back();}

    iterator begin(){return residues_.begin();}
    iterator end()  {return residues_.end();}
    const_iterator cbegin() const {return residues_.cbegin();}
    const_iterator cend()   const {return residues_.cend();}

  private:

    container_type residues_;
};

template<typename traitsT>
void PDBChain<traitsT>::push_back(const residue_type& res)
{
    if(not this->empty() && (res.chain_id() != this->chain_id()))
        throw std::invalid_argument("push invalid residue into chain");
    this->residues_.push_back(res);
    return;
}

template<typename traitsT>
void PDBChain<traitsT>::emplace_back(residue_type&& res)
{
    if(not this->empty() && (res.chain_id() != this->chain_id()))
        throw std::invalid_argument("emplace invalid residue into chain");
    this->residues_.emplace_back(std::forward<residue_type>(res));
    return;
}

}//jarngreipr
#endif /* JARNGREIPR_IO_PDB_CHAIN */
