#ifndef MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#define MJOLNIR_CORE_NAIVE_PAIR_CALCULATION_HPP
#include <mjolnir/core/SpatialPartitionBase.hpp>
#include <mjolnir/util/empty.hpp>

namespace mjolnir
{

template<typename traitsT, typename potentialT>
class NaivePairCalculation final
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

    NaivePairCalculation() = default;
    ~NaivePairCalculation() override {}
    NaivePairCalculation(NaivePairCalculation const&)            = default;
    NaivePairCalculation& operator=(NaivePairCalculation const&) = default;
    NaivePairCalculation(NaivePairCalculation&&)                 = default;
    NaivePairCalculation& operator=(NaivePairCalculation&&)      = default;

    bool valid() const noexcept override {return true;}

    void initialize(neighbor_list_type& neighbor, const system_type& sys,
                    const parameter_list_type& params) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
        MJOLNIR_LOG_FUNCTION_DEBUG();
        this->make(neighbor, sys, params);
        return;
    }

    void make(neighbor_list_type& neighbor, const system_type& sys,
              const parameter_list_type& params) override;

    bool reduce_margin(neighbor_list_type&, const real_type /*dmargin*/,
                       const system_type&, const parameter_list_type&) override
    {
        return false;
    }
    bool scale_margin(neighbor_list_type&, const real_type /*scale*/,
                      const system_type&, const parameter_list_type&) override
    {
        return false;
    }

    real_type cutoff() const noexcept override {return std::numeric_limits<real_type>::infinity();}
    real_type margin() const noexcept override {return std::numeric_limits<real_type>::infinity();}

    base_type* clone() const override {return new NaivePairCalculation();}
};

template<typename traitsT, typename potentialT>
void NaivePairCalculation<traitsT, potentialT>::make(
        neighbor_list_type& neighbors, const system_type&,
        const parameter_list_type& params)
{
    MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();
    MJOLNIR_LOG_FUNCTION_DEBUG();

    neighbors.clear();

    const auto leading_participants = params.leading_participants();
    MJOLNIR_LOG_DEBUG("leading participants size is ", leading_participants.size());

    std::vector<neighbor_type> partners;
    for(std::size_t idx=0; idx<leading_participants.size(); ++idx)
    {
        partners.clear();
        const auto i = leading_participants[idx];
        for(const auto j : params.possible_partners_of(idx, i))
        {
            if(params.has_interaction(i, j)) // likely
            {
                partners.emplace_back(j, params.prepare_params(i, j));
            }
        }
        neighbors.add_list_for(i, partners.begin(), partners.end());
    }
    return;
}

} // mjolnir
#endif /*MJOLNIR_CORE_NAIVE_PAIR_CALCULATION*/
