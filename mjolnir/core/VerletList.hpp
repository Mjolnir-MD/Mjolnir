#ifndef MJOLNIR_CORE_VERLET_LIST
#define MJOLNIR_CORE_VERLET_LIST
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <algorithm>
#include <limits>

namespace mjolnir
{

template<typename traitsT>
class VerletList
{
  public:
    typedef traitsT traits_type;
    typedef System<traits_type> system_type;
    typedef typename traits_type::boundary_type   boundary_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

    typedef ExclusionList exclusion_list_type;
    typedef NeighborList  neighbor_list_type;
    typedef neighbor_list_type::range_type range_type;

  public:

    VerletList() : margin_(0.5), current_margin_(-1.0){}
    VerletList(const real_type mgn): margin_(mgn), current_margin_(-1.0){}

    ~VerletList() = default;
    VerletList(VerletList const&) = default;
    VerletList(VerletList &&)     = default;
    VerletList& operator=(VerletList const&) = default;
    VerletList& operator=(VerletList &&)     = default;

    bool valid() const noexcept
    {
        return current_margin_ >= 0.0;
    }

    template<typename PotentialT>
    void initialize(const system_type& sys, const PotentialT& pot)
    {
        this->set_cutoff(pot.max_cutoff_length());
        this->exclusion_.make(sys, pot);
        this->make(sys);
        return;
    }

    template<typename PotentialT>
    void reconstruct(const system_type& sys, const PotentialT& pot)
    {
        this->initialize(sys, pot); // do the same thing as `initialize`
        return;
    }

    void make  (const system_type& sys);
    void update(const system_type& sys);

    real_type cutoff() const noexcept {return this->cutoff_;}
    real_type margin() const noexcept {return this->margin_;}

    range_type partners(std::size_t i) const noexcept {return neighbors_[i];}

  private:

    void set_cutoff(const real_type c) noexcept {this->cutoff_ = c;}
    void set_margin(const real_type m) noexcept {this->margin_ = m;}

  private:

    real_type      cutoff_;
    real_type      margin_;
    real_type      current_margin_;

    exclusion_list_type exclusion_;
    neighbor_list_type  neighbors_;
};

template<typename traitsT>
void VerletList<traitsT>::make(const system_type& sys)
{
    this->neighbors_.clear();

    const real_type rc = cutoff_ * (1. + margin_);
    const real_type rc2 = rc * rc;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& ri = sys[i].position;

        std::vector<std::size_t> partners;
        for(std::size_t j=i+1; j<sys.size(); ++j)
        {
            if(this->exclusion_.is_excluded(i, j))
            {
                continue;
            }

            const auto& rj = sys[j].position;
            if(length_sq(sys.adjust_direction(rj - ri)) < rc2)
            {
                partners.push_back(j);
            }
        }
        this->neighbors_.add_list_for(i, partners.begin(), partners.end());
    }
    this->current_margin_ = cutoff_ * margin_;
    return ;
}

template<typename traitsT>
void VerletList<traitsT>::update(const system_type& sys)
{
    this->current_margin_ -= sys.largest_displacement() * 2;
    if(this->current_margin_ < 0)
    {
        this->make(sys);
    }
    return ;
}

} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
