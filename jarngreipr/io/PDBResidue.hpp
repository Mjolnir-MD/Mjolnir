#ifndef JARNGREIPR_IO_PDB_RESIDUE
#define JARNGREIPR_IO_PDB_RESIDUE
#include <jarngreipr/io/PDBAtom.hpp>
#include <vector>

namespace mjolnir
{

template<typename coordT>
class PDBResidue
{
  public:
    typedef PDBAtom<coordT> atom_type;
    typedef typename atom_type::int_type   int_type;
    typedef typename atom_type::real_type  real_type;
    typedef typename atom_type::coordinate_type coordinate_type;

    typedef std::vector<atom_type> container_type;
    typedef typename container_type::value_type     value_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

    PDBResidue() = default;
    ~PDBResidue() = default;
    PDBResidue(PDBResidue const&) = default;
    PDBResidue(PDBResidue&&)      = default;
    PDBResidue& operator=(PDBResidue const&) = default;
    PDBResidue& operator=(PDBResidue&&)      = default;

    int_type           residue_id()   const noexcept
    {return this->atoms_.front().residue_id;}
    std::string const& residue_name() const noexcept
    {return this->atoms_.front().residue_name;}
    char chain_id() const noexcept
    {return this->atoms_.front().chain_id;}

    void push_back(atom_type const& atom)
    {
        if(false == is_in_same_residue(atom))
        {
            throw std::invalid_argument("mjolnir::PDBResidue::push_back:"
                "pushing an atom in different residue");
        }
        atoms_.push_back(atom);
        return;
    }
    void push_back(atom_type&& atom)
    {
        if(false == is_in_same_residue(atom))
        {
            throw std::invalid_argument("mjolnir::PDBResidue::push_back:"
                "pushing an atom in different residue");
        }
        atoms_.push_back(std::move(atom));
        return;
    }

    std::size_t size() const noexcept {return atoms_.size();}
    bool       empty() const noexcept {return atoms_.empty();}
    void clear(){atoms_.clear();}
    void resize (const std::size_t s){atoms_.resize(s);}
    void reserve(const std::size_t s){atoms_.reserve(s);}

    atom_type &      at(const std::size_t i)       {return atoms_.at(i);}
    atom_type const& at(const std::size_t i) const {return atoms_.at(i);}
    atom_type &      operator[](const std::size_t i)       noexcept
    {return atoms_[i];}
    atom_type const& operator[](const std::size_t i) const noexcept
    {return atoms_[i];}
    atom_type &      front()       noexcept {return atoms_.front();}
    atom_type const& front() const noexcept {return atoms_.front();}
    atom_type &      back()        noexcept {return atoms_.back();}
    atom_type const& back()  const noexcept {return atoms_.back();}
    iterator       begin()         noexcept {return atoms_.begin();}
    iterator       end()           noexcept {return atoms_.end();}
    const_iterator begin()   const noexcept {return atoms_.begin();}
    const_iterator end()     const noexcept {return atoms_.end();}
    const_iterator cbegin()  const noexcept {return atoms_.cbegin();}
    const_iterator cend()    const noexcept {return atoms_.cend();}

    container_type&       atoms()       noexcept {return atoms_;}
    container_type const& atoms() const noexcept {return atoms_;}

    bool is_in_same_residue(const atom_type& a) const noexcept
    {
        return atoms_.empty() ||
               (this->residue_id()   == a.residue_id &&
                this->residue_name() == a.residue_name &&
                this->chain_id()     == a.chain_id);
    }

  private:
    container_type atoms_;
};

} // jarngreipr
#endif /* JARNGREIPR_IO_PDB_RESIDUE */
