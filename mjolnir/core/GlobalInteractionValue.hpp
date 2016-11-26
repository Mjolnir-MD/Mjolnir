#ifndef MJOLNIR_GLOBAL_INTERACTION_VALUE
#define MJOLNIR_GLOBAL_INTERACTION_VALUE
#include <vector>

namespace mjolnir
{

template<typename valueT>
class GlobalInteractionValue
{
  public:
    typedef typename std::vector<std::size_t>::iterator iterator;
    typedef typename std::vector<std::size_t>::const_iterator const_iterator;

  public:
    GlobalInteractionValue() = default;
    GlobalInteractionValue(const std::size_t size): values_(size){}
    ~GlobalInteractionValue() = default;

    void emplace(const std::pair<std::size_t, valueT>& val);

    iterator begin(){return indices_.begin();}
    iterator end(){return indices_.end();}
    const_iterator cbegin() const {return indices_.cbegin();}
    const_iterator cend() const {return indices_.cend();}

    valueT&       at(const std::size_t i)       {return values_.at(i);}
    valueT const& at(const std::size_t i) const {return values_.at(i);}
    valueT&       operator[](const std::size_t i)       {return values_[i];}
    valueT const& operator[](const std::size_t i) const {return values_[i];}

  private:
    std::vector<std::size_t> indices_;
    std::vector<valueT>      values_;
};

template<typename valueT>
void GlobalInteractionValue<valueT>::emplace(
        const std::pair<std::size_t, valueT>& val)
{
    indices_.push_back(val.first);
    values_.at(val.first) = val.second;
    return ;
}

} // mjolnir
#endif /* MJOLNIR_GLOBAL_INTERACTION_VALUE */
