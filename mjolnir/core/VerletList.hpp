#ifndef MJOLNIR_CORE_VERLET_LIST
#define MJOLNIR_CORE_VERLET_LIST
#include <vector>
#include <unordered_map>
#include "SpatialPartition.hpp"
#include "BoundaryCondition.hpp"

namespace mjolnir
{

template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class VerletList : public SpatialPartition<traitsT>
{
  public:

    typedef boundaryT boundary_type;
    typedef SpatialPartition<traitsT> base_type;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::index_list index_list;

  public:

    VerletList() = default;
    VerletList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}
    VerletList(const real_type cutoff, const real_type mergin,
               const time_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}
    ~VerletList() = default;

    bool valid() const noexcept override
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void make  (const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon,
                const time_type dt) override;

    real_type cutoff() const {return this->cutoff_;}
    real_type mergin() const {return this->mergin_;}

    void set_cutoff(const real_type c){return this->cutoff_ = c;}
    void set_mergin(const real_type m){return this->mergin_ = m;}

  private:

    real_type        dt_;
    real_type        cutoff_;
    real_type        mergin_;
    real_type        current_mergin_;
};


template<typename traitsT, typename boundaryT>
void VerletList<traitsT, boundaryT>::make(const particle_container_type& pcon)
{
    this->list_.clear();
    this->list_.resize(pcon.size());

    const real_type rc = (cutoff_ + mergin_);
    const real_type rc2 = rc * rc;
    for(std::size_t i=0; i<pcon.size()-1; ++i)
    {
        const coordinate_type& ri = pcon[i].position;

        const auto cbeg = this->except_.at(i).cbegin();
        const auto cend = this->except_.at(i).cend();
        for(std::size_t j=i+1; j<pcon.size(); ++j)
        {
            if((std::find(cbeg, cend, j) == cend) &&
               (length_sq(boundary_type::adjust_direction(pcon[j].position - ri))
                < rc2))
                this->list_.at(i).push_back(j);
        }
    }
    this->current_mergin_ = mergin_;
    return ;
}

template<typename traitsT, typename boundaryT>
void VerletList<traitsT, boundaryT>::update(const particle_container_type& pcon)
{
    real_type max_speed = 0.;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
        max_speed = std::max(max_speed, length_sq(iter->velocity));

    this->current_mergin_ -= std::sqrt(max_speed) * dt_ * 2.;
    if(this->current_mergin_ < 0.)
        this->make(pcon);

    return ;
}


template<typename traitsT, typename boundaryT>
void VerletList<traitsT, boundaryT>::update(const particle_container_type& pcon,
        const time_type dt)
{
    this->dt_ = dt;
    this->update(pcon);
    return ;
}


} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
