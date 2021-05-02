#ifndef MJOLNIR_CORE_VERLET_LIST_HPP
#define MJOLNIR_CORE_VERLET_LIST_HPP
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <algorithm>
#include <limits>

namespace mjolnir
{

template<typename traitsT, typename potentialT>
class VerletList final
    : public SpatialPartitionBase<traitsT, potentialT>
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using base_type           = SpatialPartitionBase<traits_type, potential_type>;

    using system_type         = typename base_type::system_type;
    using boundary_type       = typename base_type::boundary_type;
    using real_type           = typename base_type::real_type;
    using coordinate_type     = typename base_type::coordinate_type;
    using neighbor_list_type  = typename base_type::neighbor_list_type;
    using neighbor_type       = typename base_type::neighbor_type;
    using range_type          = typename base_type::range_type;
    using parameter_list_type = typename base_type::parameter_list_type;

  public:

    VerletList() : margin_(0.5), current_margin_(-1.0){}
    explicit VerletList(const real_type mgn): margin_(mgn), current_margin_(-1.0){}

    ~VerletList() override {}
    VerletList(VerletList const&) = default;
    VerletList(VerletList &&)     = default;
    VerletList& operator=(VerletList const&) = default;
    VerletList& operator=(VerletList &&)     = default;

    bool valid() const noexcept override
    {
        return current_margin_ >= 0.0;
    }

    void initialize(neighbor_list_type& neighbors, const system_type& sys,
                    const parameter_list_type& params) override
    {
        this->set_cutoff(params.max_cutoff_length());
        this->make(neighbors, sys, params);
        return;
    }

    void make(neighbor_list_type& neighbors, const system_type& sys,
              const parameter_list_type& params) override;

    bool reduce_margin(neighbor_list_type& neighbors, const real_type dmargin,
                       const system_type& sys, const parameter_list_type& params) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, params);
            return true;
        }
        return false;
    }
    bool scale_margin(neighbor_list_type& neighbors, const real_type scale,
                      const system_type& sys, const parameter_list_type& params) override
    {
        this->current_margin_ = (cutoff_ + current_margin_) * scale - cutoff_;
        if(this->current_margin_ < 0)
        {
            this->make(neighbors, sys, params);
            return true;
        }
        return false;
    }

    real_type cutoff() const noexcept override {return this->cutoff_;}
    real_type margin() const noexcept override {return this->margin_;}

    base_type* clone() const override
    {
        return new VerletList(margin_);
    }

  private:

    void set_cutoff(const real_type c) noexcept {this->cutoff_ = c;}
    void set_margin(const real_type m) noexcept {this->margin_ = m;}

  private:

    real_type cutoff_;
    real_type margin_;
    real_type current_margin_;
};

template<typename traitsT, typename potentialT>
void VerletList<traitsT, potentialT>::make(
        neighbor_list_type& neighbors, const system_type& sys,
        const parameter_list_type& params)
{
    neighbors.clear();

    const real_type rc = cutoff_ * (1. + margin_);
    const real_type rc2 = rc * rc;

    const auto leading_participants = params.leading_participants();

    std::vector<neighbor_type> partner;
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        partner.clear();
        const auto   i = leading_participants[idx];
        const auto& ri = sys.position(i);

        for(const auto j : params.possible_partners_of(idx, i))
        {
            if(!params.has_interaction(i, j))
            {
                continue;
            }
            const auto& rj = sys.position(j);
            if(math::length_sq(sys.adjust_direction(ri, rj)) < rc2)
            {
                partner.emplace_back(j, potential_type(params.prepare_params(i, j)));
            }
        }
        // because j is searched sequencially, sorting is not needed.
        neighbors.add_list_for(i, partner.begin(), partner.end());
    }
    this->current_margin_ = cutoff_ * margin_;
    return ;
}
} // mjolnir
#endif/* MJOLNIR_CORE_VERLET_LIST */
