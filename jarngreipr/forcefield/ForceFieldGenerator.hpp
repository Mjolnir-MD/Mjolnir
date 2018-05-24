#ifndef JARNGREIPR_FORCEFIELD_GENERATOR
#define JARNGREIPR_FORCEFIELD_GENERATOR
#include <jarngreipr/model/Bead.hpp>
#include <extlib/toml/toml.hpp>

namespace jarngreipr
{

template<typename realT>
class ForceFieldGenerator
{
  public:
    typedef realT real_type;
    typedef Bead<real_type> bead_type;

  public:
    virtual ~IntraChainForceFieldGenerator() = default;

    //!@brief generate forcefield parameter values
    virtual void generate(toml::Table& out,
            const std::vector<std::vector<std::unique_ptr<bead_type>>>& chains
            ) const = 0;

    //!@brief if chain contains invalid bead, return false.
    virtual bool check_beads_kind(
            const std::vector<std::unique_ptr<bead_type>>& chain) const = 0;
};

} // mjolnir
#endif// JARNGREIPR_FORCEFIELD_GENERATOR
