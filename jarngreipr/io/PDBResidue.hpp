#ifndef JARNGREIPR_IO_PDB_RESIDUE
#define JARNGREIPR_IO_PDB_RESIDUE
#include "PDBAtom.hpp"
#include <vector>
#include <stdexcept>

namespace jarngreipr
{

template<typename traitsT>
class PDBResidue
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef PDBAtom<traits_type> atom_type;
    typedef typename atom_type::int_type int_type;
    typedef std::vector<atom_type> container_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

    PDBResidue() = default;
    ~PDBResidue() = default;
    PDBResidue(PDBResidue const&) = default;
    PDBResidue(PDBResidue&&)      = default;
    PDBResidue& operator=(PDBResidue const&) = default;
    PDBResidue& operator=(PDBResidue&&)      = default;

    int_type           residue_id()   const {return atoms_.front().residue_id;}
    std::string const& residue_name() const {return atoms_.front().residue_name;}
    std::string const& chain_id()     const {return atoms_.front().chain_id;}

    void push_back(const atom_type& atom);
    void emplace_back(atom_type&& atom);

    std::size_t size() const {return atoms_.size();}
    bool       empty() const {return atoms_.empty();}
    void clear(){atoms_.clear();}
    void resize(const std::size_t s){atoms_.reserve(s);}
    void reserve(const std::size_t s){atoms_.reserve(s);}

    atom_type &      at(const std::size_t i)       {return atoms_.at(i);}
    atom_type const& at(const std::size_t i) const {return atoms_.at(i);}
    atom_type &      operator[](const std::size_t i)       {return atoms_[i];}
    atom_type const& operator[](const std::size_t i) const {return atoms_[i];}
    atom_type &      front()       {return atoms_.front();}
    atom_type const& front() const {return atoms_.front();}
    atom_type &      back()        {return atoms_.back();}
    atom_type const& back()  const {return atoms_.back();}
    iterator begin() {return atoms_.begin();}
    iterator end()   {return atoms_.end();}
    const_iterator cbegin() const {return atoms_.cbegin();}
    const_iterator cend()   const {return atoms_.cend();}

  private:
    container_type atoms_;
};

template<typename traitsT>
void PDBResidue<traitsT>::push_back(const atom_type& atom)
{
    if(not atoms_.empty() && (this->residue_id() != atom.residue_id ||
                              this->residue_name() != atom.residue_name ||
                              this->chain_id() != atom.chain_id))
        throw std::invalid_argument("push invalid atom to residue");
    atoms_.push_back(atom);
    return;
}

template<typename traitsT>
void PDBResidue<traitsT>::emplace_back(atom_type&& atom)
{
    if(not atoms_.empty() && (this->residue_id() != atom.residue_id ||
                              this->residue_name() != atom.residue_name ||
                              this->chain_id() != atom.chain_id))
        throw std::invalid_argument("emplace invalid atom to residue");
    atoms_.emplace_back(std::forward<atom_type>(atom));
    return;
}

} // jarngreipr
#endif /* JARNGREIPR_IO_PDB_RESIDUE */
