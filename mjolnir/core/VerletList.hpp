#ifndef MJOLNIR_CORE_VERLET_LIST
#define MJOLNIR_CORE_VERLET_LIST
#include <vector>
#include <unordered_map>
#include "SpatialPartition.hpp"

namespace mjolnir
{

template<typename traitsT>
class VerletList : public SpatialPartition<traitsT>
{
  public:

    typedef SpatialPartition<traitsT> base_type;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::index_list index_list;
    typedef std::vector<index_list> verlet_list_type;

  public:

    VerletList() = default;
    VerletList(const real_type cutoff, const real_type mergin, const time_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}
    ~VerletList() = default;

    bool valid() const noexcept override {return current_mergin_ >= 0.;}

    void make  (const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon, const time_type dt) override;

    real_type const& cutoff() const {return this->cutoff_;}
    real_type &      cutoff()       {return this->cutoff_;}
    real_type const& mergin() const {return this->mergin_;}
    real_type &      mergin()       {return this->mergin_;}

    index_list const& partners(const std::size_t i) const override {return list_.at(i);}
    index_list &      partners(const std::size_t i)       override {return list_.at(i);}

  private:

    real_type        dt_;
    real_type        cutoff_;
    real_type        mergin_;
    real_type        current_mergin_;
    verlet_list_type list_;
};

template<typename traitsT>
void VerletList<traitsT>::make(const particle_container_type& pcon)
{
    list_.clear();
    list_.resize(pcon.size());

    const real_type rc = (cutoff_ + mergin_);
    const real_type rc2 = rc * rc;
    for(std::size_t i=0; i<pcon.size()-1; ++i)
    {
        const coordinate_type ri = pcon[i].position;

        for(std::size_t j=i+1; j<pcon.size(); ++j)
        {
            const coordinate_type rij = pcon[j].position - ri;
            const real_type r2 = length_sq(rij);
            if(r2 < rc2)
            {
                list_.at(i).push_back(j);
            }
        }
    }
    this->current_mergin_ = mergin_;
    return ;
}

template<typename traitsT>
void VerletList<traitsT>::update(const particle_container_type& pcon)
{
    real_type max_speed = 0.;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
        max_speed = std::max(max_speed, length_sq(iter->velocity));

    this->current_mergin_ -= std::sqrt(max_speed) * dt_ * 2.;
    if(this->current_mergin_ < 0.)
        this->make(pcon);

    return ;
}


template<typename traitsT>
void VerletList<traitsT>::update(const particle_container_type& pcon,
        const time_type dt)
{
    real_type max_speed = 0.;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
        max_speed = std::max(max_speed, length_sq(iter->velocity));

    const real_type reduction = std::sqrt(max_speed) * dt * 2;

    this->current_mergin_ -= reduction;
    if(this->current_mergin_ < 0.)
        this->make(pcon);

    return ;
}


} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
