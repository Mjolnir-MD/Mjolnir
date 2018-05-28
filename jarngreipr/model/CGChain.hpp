#ifndef JARNGREIPR_MODEL_CG_CHAIN_H
#define JARNGREIPR_MODEL_CG_CHAIN_H
#include <jarngreipr/memory/Beads.hpp>
#include <vector>
#include <memory>

namespace jarngreipr
{

template<typename realT>
class CGChain
{
  public:
    typedef realT real_type;
    typedef Bead<real_type> bead_type;
    typedef std::shared_ptr<bead_type> bead_ptr;
    typedef std::vector<bead_ptr> container_type;

  public:

    CGChain(const std::string& name): name_(name){}

    CGChain(const CGChain&) = default;
    CGChain(CGChain&&)      = default;
    CGChain& operator=(const CGChain&) = default;
    CGChain& operator=(CGChain&&)      = default;
    ~CGChain() = default;

    container_type const& beads() const noexcept {return beads_;}
    container_type&       beads()       noexcept {return beads_;}
    std::string const&    name()  const noexcept {return name_;}
    std::string&          name()        noexcept {return name_;}

  protected:

    std::string    name_;
    container_type beads_;
};

}//jarngreipr
#endif // JARNGREIPR_MODEL_CG_CHAIN_H
