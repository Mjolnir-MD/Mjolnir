#ifndef JARNGREIPR_BEAD
#define JARNGREIPR_BEAD
#include <jarngreipr/io/PDBAtom.hpp>
#include <vector>
#include <string>

namespace jarngreipr
{

template<typename traitsT>
class Bead
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef PDBAtom<traits_type> atom_type;
    typedef std::vector<atom_type> container_type;

  public:

    Beads() = default;
    explicit Beads(const container_type& atoms) : atoms_(atoms){}
    explicit Beads(const std::string& name) : name_(name){}
    Beads(const container_type& atoms, const std::string& name)
        : atoms_(atoms), name_(name){}
    ~CarbonAlpha() = default;

    virtual coordinate_type position(const std::size_t i) const = 0;

    residue_type const& atoms() const {return atoms_;}
    residue_type &      atoms()       {return atoms_;}

    std::string const& name() const {return name_;}
    std::string &      name()       {return name_;}

  protected:

    std::string     name_;
    container_type  atoms_;
};

}//jarngreipr
#endif /* JARNGREIPR_BEAD */
