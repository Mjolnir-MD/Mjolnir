#ifndef JARNGREIPR_MODEL_BEAD
#define JARNGREIPR_MODEL_BEAD
#include <jarngreipr/io/PDBAtom.hpp>
#include <vector>
#include <string>

namespace mjolnir
{

template<typename coordT>
class Bead
{
  public:
    typedef coordT                 coord_type;
    typedef PDBAtom<coord_type>    atom_type;
    typedef std::vector<atom_type> container_type;

  public:
    Bead() = default;
    virtual ~Bead() = default;

    explicit Bead(const container_type& atoms) : atoms_(atoms){}
    explicit Bead(container_type&&      atoms) : atoms_(std::move(atoms)){}
    explicit Bead(const std::string& name)     : name_(name){}
    explicit Bead(std::string&& name)          : name_(std::move(name)){}

    Bead(const container_type& atoms, const std::string& name)
        : atoms_(atoms), name_(name)
    {}
    Bead(const container_type& atoms, std::string&& name)
        : atoms_(atoms), name_(std::move(name))
    {}
    Bead(container_type&& atoms, const std::string& name)
        : atoms_(std::move(atoms)), name_(name)
    {}
    Bead(container_type&& atoms, std::string&& name)
        : atoms_(std::move(atoms)), name_(std::move(name))
    {}

    virtual coord_type  position() const = 0;
    virtual std::string attribute(const std::string& attr_name) const = 0;

    void assign(const atom_type& atom) {atoms_.push_back(atom);}

    container_type const& atoms() const noexcept {return atoms_;}
    container_type &      atoms()       noexcept {return atoms_;}

    std::string const& name() const noexcept {return name_;}
    std::string &      name()       noexcept {return name_;}

  protected:

    std::string     name_;
    container_type  atoms_;
};

}//jarngreipr
#endif /* JARNGREIPR_BEAD */
