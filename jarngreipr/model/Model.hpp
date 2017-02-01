#ifndef JARNGREIPR_MODEL
#define JARNGREIPR_MODEL
#include "Bead.hpp"

namespace jarngreipr
{

template<typename traitsT>
class Model
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef PDBChain<traits_type> chain_type;
    typedef std::vector<chain_type> chain_container_type;
    typedef Bead<traits_type> bead_type;
    typedef std::vector<bead_type> bead_container_type;

    template<std::size_t N>
    using interaction_type = std::pair<std::array<std::size_t, N>, real_type>;

  public:

    Model() = default;
    ~Model() = default;

    virtual void make(const chain_container_type& chain) = 0;

    bead_container_type&       beads()       {return beads_;}
    bead_container_type const& beads() const {return beads_;}

  protected:

    bead_container_type beads_;
};

}//jarngreipr
#endif /* JARNGREIPR_MODEL */
