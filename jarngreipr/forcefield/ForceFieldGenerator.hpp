#ifndef JARNGREIPR_FORCEFIELD_GENERATOR
#define JARNGREIPR_FORCEFIELD_GENERATOR
#include <jarngreipr/model/CGChain.hpp>
#include <extlib/toml/toml.hpp>
#include <memory>

namespace jarngreipr
{

template<typename realT>
class ForceFieldGenerator
{
  public:
    typedef realT real_type;
    typedef Bead<real_type>    bead_type;
    typedef CGChain<real_type> chain_type;

  public:
    virtual ~ForceFieldGenerator() = default;

    //!@brief generate forcefield parameter values
    virtual void generate(toml::Table& out,
        const std::vector<chain_type>& chains) const = 0;

    //!@brief generate inter-chain parameters if it's defined.
    virtual void generate(toml::Table& out,
        const std::vector<chain_type>& lhs,
        const std::vector<chain_type>& rhs) const = 0;

    //!@brief if chain contains invalid bead, return false.
    virtual bool check_beads_kind(const chain_type& chain) const = 0;
};

} // mjolnir
#endif// JARNGREIPR_FORCEFIELD_GENERATOR
