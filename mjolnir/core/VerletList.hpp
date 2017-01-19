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
    typedef std::vector<index_list> except_list_type;

  public:

    VerletList() = default;
    VerletList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}
    VerletList(const real_type cutoff, const real_type mergin, const time_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.)
    {}
    ~VerletList() = default;

    bool valid() const noexcept override {return current_mergin_ >= 0. || dt_ == 0.;}

    void make  (const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon, const time_type dt) override;

    real_type const& cutoff() const {return this->cutoff_;}
    real_type &      cutoff()       {return this->cutoff_;}
    real_type const& mergin() const {return this->mergin_;}
    real_type &      mergin()       {return this->mergin_;}

    void add_except(const index_type i, const index_type j);
    void set_except(const except_list_type& ex){except_ = ex;}
    void set_except(except_list_type&& ex){except_ = std::forward<except_list_type>(ex);}

    index_list const& partners(const std::size_t i) const override {return list_[i];}
    index_list &      partners(const std::size_t i)       override {return list_[i];}

  private:

    real_type        dt_;
    real_type        cutoff_;
    real_type        mergin_;
    real_type        current_mergin_;
    verlet_list_type list_;
    except_list_type except_;
};

template<typename traitsT>
inline void
VerletList<traitsT>::add_except(const index_type i, const index_type j)
{
    const index_type acc = std::min(i, j);
    const index_type dnr = std::max(i, j);

    if(this->except_.size() < acc)
        this->except_.resize(acc);

    this->except_.at(acc).push_back(dnr);
    return ;
}


template<typename traitsT>
void VerletList<traitsT>::make(const particle_container_type& pcon)
{
    list_.clear();
    list_.resize(pcon.size());

    const real_type rc = (cutoff_ + mergin_);
    const real_type rc2 = rc * rc;
    if(except_.empty())
    {
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
    }
    else
    {
        for(std::size_t i=0; i<pcon.size()-1; ++i)
        {
            const coordinate_type ri = pcon[i].position;

            if(except_.at(i).empty())
            {
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
            else
            {
                const auto cbeg = except_.at(i).cbegin();
                const auto cend = except_.at(i).cend();
                for(std::size_t j=i+1; j<pcon.size(); ++j)
                {
                    if(std::find(cbeg, cend, j) != cend)
                        continue;
                    const coordinate_type rij = pcon[j].position - ri;
                    const real_type r2 = length_sq(rij);
                    if(r2 < rc2)
                    {
                        list_.at(i).push_back(j);
                    } 
                }
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
    this->dt_ = dt;
    this->update(pcon);
    return ;
}


} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
