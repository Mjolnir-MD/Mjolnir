#ifndef MJOLNIR_CORE_NAIVE_PAIR_CALCULATION
#define MJOLNIR_CORE_NAIVE_PAIR_CALCULATION
#include "SpatialPartition.hpp"
#include <iostream>

namespace mjolnir
{

template<typename traitsT>
class NaivePairCalculation : public SpatialPartition<traitsT>
{
  public:

    typedef SpatialPartition<traitsT> base_type;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::index_list index_list;
    typedef std::vector<index_list> list_type;

  public:

    NaivePairCalculation(): dirty_(true){};
    ~NaivePairCalculation() = default;

    NaivePairCalculation(NaivePairCalculation&&) = default;
    NaivePairCalculation(const NaivePairCalculation&) = delete;
    NaivePairCalculation& operator=(NaivePairCalculation&&) = default;
    NaivePairCalculation& operator=(const NaivePairCalculation&) = delete;

    bool valid() const noexcept override {return true;}

    void make  (const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon,
                const time_type dt) override;

    void add_except(const index_type i, const index_type j);
    void set_except(const list_type& ex){excepts_ = ex;}
    void set_except(list_type&& ex){excepts_ = std::forward<list_type>(ex);}

    index_list const& partners(const std::size_t i) const override {return list_.at(i);}
    index_list &      partners(const std::size_t i)       override {return list_.at(i);}

  private:

    bool dirty_;
    list_type list_;
    list_type excepts_;
};

template<typename traitsT>
inline void NaivePairCalculation<traitsT>::add_except(
        const index_type i, const index_type j)
{
    const index_type acc = std::min(i, j);
    const index_type dnr = std::max(i, j);

    if(this->excepts_.size() < acc)
        this->excepts_.resize(acc);

    this->excepts_.at(acc).push_back(dnr);
    return ;
}

template<typename traitsT>
void NaivePairCalculation<traitsT>::make(const particle_container_type& pcon)
{
    list_.resize(pcon.size());
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        if(excepts_.at(i).empty())
        {
            list_.at(i).resize(pcon.size() - i - 1);
            for(std::size_t j=i+1; j<pcon.size(); ++j)
                list_.at(i).at(j-i+1) = j;
        }
        else
        {
            list_.at(i).reserve(pcon.size() - i - 1 - excepts_.at(i).size());
            const auto cbeg = excepts_.at(i).cbegin();
            const auto cend = excepts_.at(i).cend();
            for(std::size_t j=i+1; j<pcon.size(); ++j)
            {
                if(std::find(cbeg, cend, j) == cend)
                    list_.at(i).push_back(j);
            }
        }
    }

    dirty_ = false;
    return;
}

template<typename traitsT>
void NaivePairCalculation<traitsT>::update(const particle_container_type& pcon)
{
    if(this->dirty_)
        this->make(pcon);
    return ;
}

template<typename traitsT>
void NaivePairCalculation<traitsT>::update(const particle_container_type& pcon,
        const time_type dt)
{
    if(this->dirty_)
        this->make(pcon);
    return ;
}

} // mjolnir
#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
