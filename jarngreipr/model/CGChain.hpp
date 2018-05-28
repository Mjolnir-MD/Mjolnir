#ifndef JARNGREIPR_MODEL_CG_CHAIN_H
#define JARNGREIPR_MODEL_CG_CHAIN_H
#include <jarngreipr/model/Bead.hpp>
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
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    CGChain(const std::string& name): name_(name){}

    CGChain(const CGChain&) = default;
    CGChain(CGChain&&)      = default;
    CGChain& operator=(const CGChain&) = default;
    CGChain& operator=(CGChain&&)      = default;
    ~CGChain() = default;

    std::size_t size() const noexcept {return beads_.size();}

    bead_ptr&       operator[](const std::size_t i)       noexcept {return beads_[i];}
    bead_ptr const& operator[](const std::size_t i) const noexcept {return beads_[i];}
    bead_ptr&       at(const std::size_t i)       {return beads_.at(i);}
    bead_ptr const& at(const std::size_t i) const {return beads_.at(i);}

    iterator       begin()        noexcept {return beads_.begin();}
    iterator       end()          noexcept {return beads_.end();}
    const_iterator begin()  const noexcept {return beads_.begin();}
    const_iterator end()    const noexcept {return beads_.end();}
    const_iterator cbegin() const noexcept {return beads_.cbegin();}
    const_iterator cend()   const noexcept {return beads_.cend();}

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
