#ifndef JARNGREIPR_FORCEFIELD_GENERATOR
#define JARNGREIPR_FORCEFIELD_GENERATOR
#include <jarngreipr/model/Bead.hpp>
#include <extlib/toml/toml.hpp>

namespace mjolnir
{

template<typename realT>
class IntraChainForceFieldGenerator
{
  public:
    typedef realT real_type;
    typedef Bead<real_type> bead_type;

  public:
    virtual ~IntraChainForceFieldGenerator() = default;

    //!@brief generate parameter values and write out to ostream.
    //@return connection_information. for scaling or ignoring Global interaction
    virtual void generate(toml::Table& out,
            const std::vector<std::unique_ptr<bead_type>>& chain) const = 0;

    //!@brief if chain contains invalid bead, return false.
    virtual bool check_beads_kind(
            const std::vector<std::unique_ptr<bead_type>>& chain) const = 0;
};

template<typename realT>
class InterChainForceFieldGenerator
{
  public:
    typedef realT real_type;
    typedef Bead<real_type> bead_type;
    typedef std::unique_ptr<bead_type> baed_ptr;

  public:
    virtual ~InterChainForceFieldGenerator() = default;

    //!@brief generate parameter values and write out to ostream.
    //@return connection_information. for scaling or ignoring Global interaction
    virtual void generate(toml::Table& out,
        const std::vector<std::vector<bead_ptr>>& chain) const = 0;

    //!@brief if chain contains invalid bead, return false.
    virtual bool check_beads_kind(
        const std::vector<std::vector<bead_type>>& chain) const = 0;
};

} // mjolnir
#endif// JARNGREIPR_FORCEFIELD_GENERATOR
