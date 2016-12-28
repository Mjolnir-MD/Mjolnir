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

    index_list const& partners(const std::size_t i) const override {return list_.at(i);}
    index_list &      partners(const std::size_t i)       override {return list_.at(i);}

  private:

    bool dirty_;
    list_type list_;
};

template<typename traitsT>
void NaivePairCalculation<traitsT>::make(const particle_container_type& pcon)
{
    list_.resize(pcon.size());
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        list_.at(i).resize(pcon.size() - i - 1);
        for(std::size_t j=i+1; j<pcon.size(); ++j)
            list_.at(i).at(j-i-1) = j;
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
