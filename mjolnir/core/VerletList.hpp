#ifndef MJOLNIR_CORE_VERLET_LIST_HPP
#define MJOLNIR_CORE_VERLET_LIST_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/NeighborList.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <algorithm>
#include <limits>

namespace mjolnir
{

template<typename traitsT, typename parameterT>
class VerletList
{
  public:
    using traits_type         = traitsT;
    using system_type         = System<traits_type>;
    using boundary_type       = typename traits_type::boundary_type;
    using real_type           = typename traits_type::real_type;
    using coordinate_type     = typename traits_type::coordinate_type;
    using exclusion_list_type = ExclusionList;

    using parameter_type      = parameterT;
    using neighbor_list_type  = NeighborList<parameter_type>;
    using neighbor_type       = typename neighbor_list_type::neighbor_type;
    using range_type          = typename neighbor_list_type::range_type;

  public:
    VerletList() : margin_(0.5), current_margin_(-1.0){}
    explicit VerletList(const real_type mgn): margin_(mgn), current_margin_(-1.0){}

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
        this->make(sys, pot);
        return;
    }

    template<typename PotentialT>
    void make  (const system_type& sys, const PotentialT& pot);

    template<typename PotentialT>
    void update(const real_type, const system_type&, const PotentialT&);

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

template<typename traitsT, typename parameterT>
template<typename potentialT>
void VerletList<traitsT, parameterT>::update(
        const real_type dmargin, const system_type& sys, const potentialT& pot)
{
    this->current_margin_ -= dmargin;
    if(this->current_margin_ < 0)
    {
        this->make(sys, pot);
    }
    return ;
}

template<typename traitsT, typename parameterT>
template<typename potentialT>
void VerletList<traitsT, parameterT>::make(
        const system_type& sys, const potentialT& pot)
{
    static_assert(std::is_same<typename potentialT::parameter_type,
        parameter_type>::value, "VerletList: invalid template argumnet: "
        "potentialT::parameter_type should be equal to verletlist::parameterT");

    this->neighbors_.clear();

    // `participants` is a list that contains indices of particles that are
    // related to the potential.
    const auto& participants = pot.participants();

    // To construct NeighborList correctly, we need to add partners to *all*
    // the particles regardless of whether the particle is a participant or not.
    //
    // So the following code first check whether an i-th particle is a
    // participant to skip needless calculation. But searching indices takes a
    // time.
    //
    // Current implementation code assumes that `participants` is sorted
    // (it IS sorted in potential::initialize()). With this assumption, we can
    // skip searching by...
    // - check the index is the same as the first element of `participant`.
    // - if they are the same, pop the first element and calc neighbors.
    // - otherwise, skip neighbor calculation and continue the loop.
    //
    // To emulate pop_front (we can't use it here because `participants` is
    // borrowed from `potential`), it uses `participant_index` to represent the
    // current `first element`.

    const real_type rc = cutoff_ * (1. + margin_);
    const real_type rc2 = rc * rc;

    std::vector<neighbor_type> partner;
    std::size_t participant_index = 0;
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        partner.clear();
        if(participant_index < participants.size() &&
           participants[participant_index] == i)
        {
            ++participant_index;
            const auto& ri = sys.position(i);

            for(std::size_t j=i+1; j<sys.size(); ++j)
            {
                if(this->exclusion_.is_excluded(i, j))
                {
                    continue;
                }

                const auto& rj = sys.position(j);
                if(math::length_sq(sys.adjust_direction(rj - ri)) < rc2)
                {
                    partner.emplace_back(j, pot.prepare_params(i, j));
                }
            }
            // because j is searched sequencially, sorting is not needed.
        }
        this->neighbors_.add_list_for(i, partner.begin(), partner.end());
    }
    this->current_margin_ = cutoff_ * margin_;
    return ;
}

} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
