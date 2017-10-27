#ifndef JARNGREIPR_MODEL_CG_CHAIN
#define JARNGREIPR_MODEL_CG_CHAIN
#include <jarngreipr/model/Bead.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <utility>
#include <memory>
#include <vector>
#include <string>

namespace mjolnir
{

template<typename coordT>
class CGChain
{
  public:
    typedef coordT       coord_type;
    typedef Bead<coordT> bead_type;
    typedef std::unique_ptr<bead_type> bead_ptr_type;
    typedef std::vector<bead_ptr_type> container_type;
    typedef typename container_type::value_type      value_type;
    typedef typename container_type::size_type       size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::allocator_type  allocator_type;
    typedef typename container_type::pointer         pointer;
    typedef typename container_type::reference       reference;
    typedef typename container_type::const_pointer   const_pointer;
    typedef typename container_type::const_reference const_reference;
    typedef typename container_type::iterator        iterator;
    typedef typename container_type::const_iterator  const_iterator;
    typedef typename container_type::reverse_iterator       reverse_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;

  public:
    CGChain()  = default;
    ~CGChain() = default;
    CGChain(CGChain const&) = default;
    CGChain(CGChain&&)      = default;
    CGChain& operator=(CGChain const&) = default;
    CGChain& operator=(CGChain&&)      = default;

    void push_back(bead_ptr_type&& bead)
    {
        this->container_.push_back(std::move(bead));
    }
    template<typename ...Ts>
    void emplace_back(Ts&& ... values)
    {
        this->container_.push_back(
                make_unique<bead_type>(std::forward<Ts>(values)...));
    }

    bool        empty() const noexcept {return container_.empty();}
    std::size_t size()  const noexcept {return container_.size();}
    void clear(){return container_.clear();}
    void reserve(const std::size_t s){container_.reserve(s);}

    bead_ptr_type&       front()       noexcept {return container_.front();}
    bead_ptr_type const& front() const noexcept {return container_.front();}
    bead_ptr_type&       back()        noexcept {return container_.back();}
    bead_ptr_type const& back()  const noexcept {return container_.back();}
    bead_ptr_type&       operator[](const std::size_t i)       noexcept {return container_[i];}
    bead_ptr_type const& operator[](const std::size_t i) const noexcept {return container_[i];}
    bead_ptr_type&       at(const std::size_t i)       {return container_.at(i);}
    bead_ptr_type const& at(const std::size_t i) const {return container_.at(i);}

    iterator       begin()        noexcept {container_.begin();}
    iterator       end()          noexcept {container_.end();}
    const_iterator begin()  const noexcept {container_.begin();}
    const_iterator end()    const noexcept {container_.end();}
    const_iterator cbegin() const noexcept {container_.cbegin();}
    const_iterator cend()   const noexcept {container_.cend();}

  private:
    container_type container_;
};

} // mjolnir
#endif // JARNGREIPR_MODEL_CG_CHAIN
