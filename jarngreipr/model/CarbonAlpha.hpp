#ifndef JARNGREIPR_MODEL_CARBON_ALPHA
#define JARNGREIPR_MODEL_CARBON_ALPHA
#include <jarngreipr/model/Bead.hpp>
#include <jarngreipr/io/PDBResidue.hpp>
#include <jarngreipr/util/string.hpp>
#include <algorithm>
#include <stdexcept>
#include <string>

namespace mjolnir
{

/*! @brief carbon alpha 1 beads per amino acid model */
template<typename coordT>
class CarbonAlpha final : public Bead<coordT>
{
  public:

    typedef Bead<coordT> base_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::atom_type       atom_type;
    typedef typename base_type::container_type  container_type;

    typedef PDBResidue<coordT> residue_type;

  public:

    CarbonAlpha() = default;
    ~CarbonAlpha() override = default;

    explicit CarbonAlpha(const residue_type& residue)
        : base_type(residue.atoms(), residue.residue_name())
    {}

    explicit CarbonAlpha(const container_type& atoms): base_type(atoms){}
    explicit CarbonAlpha(container_type&& atoms) : base_type(std::move(atoms)){}
    explicit CarbonAlpha(const std::string& name): base_type(name){}
    explicit CarbonAlpha(std::string&& name)     : base_type(std::move(name)){}

    CarbonAlpha(const residue_type& residue, const std::string& name)
        : base_type(residue.atoms(), name)
    {}
    CarbonAlpha(const container_type&& atoms, const std::string& name)
        : base_type(atoms, name)
    {}
    CarbonAlpha(container_type&& atoms, const std::string& name)
        : base_type(std::move(atoms), name)
    {}
    CarbonAlpha(const container_type& atoms, std::string&& name)
        : base_type(atoms, std::move(name))
    {}
    CarbonAlpha(container_type&& atoms, std::string&& name)
        : base_type(std::move(atoms), std::move(name))
    {}

    coordinate_type position() const override;

    std::string     attribute(const std::string& n) const override
    {return ""_str;}

    std::string kind() const override {return "CarbonAlpha"_str;}
};

template<typename coordT>
typename CarbonAlpha<coordT>::coordinate_type
CarbonAlpha<coordT>::position() const
{
    const auto finder = [](const atom_type& a){return a.atom_name == "CA"_str;};
    const std::size_t num_ca =
        std::count_if(this->atoms_.cbegin(), this->atoms_.cend(), finder);
    if(num_ca == 0)
    {
        throw std::runtime_error("mjolnir::model::CarbonAlpha::position: "
                "no c-alpha atom in this residue");
    }
    else if(num_ca != 1)
    {
        throw std::runtime_error("mjolnir::model::CarbonAlpha::position: "_str +
                "multiple ("_str + std::to_string(num_ca) +
                ") C-alphas exist in this beads "_str + this->name_);
    }
    const auto ca =
        std::find_if(this->atoms_.cbegin(), this->atoms_.cend(), finder);
    return ca->position;
}

}//jarngreipr
#endif /*JARNGREIPR_CARBON_ALPHA*/
