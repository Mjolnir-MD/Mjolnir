#ifndef JARNGREIPR_CARBON_ALPHA
#define JARNGREIPR_CARBON_ALPHA
#include "Bead.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

namespace jarngreipr
{

/*! @brief carbon alpha 1 beads model */
template<typename traitsT>
class CarbonAlpha : public Bead<traitsT>
{
  public:

    typedef Bead<traitsT> base_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::container_type container_type;
    typedef typename base_type::atom_type atom_type;

  public:

    CarbonAlpha() = default;
    CarbonAlpha(const residue_type& residue) : base_type(residue.atoms()){}
    CarbonAlpha(const container_type& atoms) : base_type(atoms){}
    ~CarbonAlpha() = default;

    coordinate_type position(const std::size_t i = 0) const override;
};

template<typename traitsT>
typename CarbonAlpha<traitsT>::coordinate_type
CarbonAlpha<traitsT>::position(const std::size_t i) const
{
    auto finder = [](const atom_type& a){return a.atom_name == "CA";};

    const std::size_t num_ca =
        std::count_if(atoms_.cbegin(), atoms_.cend(), finder);
    if(num_ca == 0)
        throw std::runtime_error("Ca: no c-alpha atom in this residue");
    else if(num_ca <= i)
        throw std::out_of_range(std::string("Ca: ") + std::to_string(i) +
                                std::string("-th Calpha does not exists"));

    const auto ca = std::find_if(atoms_.cbegin(), atoms_.cend(), finder);
    return ca->position;
}

}//jarngreipr
#endif /*JARNGREIPR_CARBON_ALPHA*/
