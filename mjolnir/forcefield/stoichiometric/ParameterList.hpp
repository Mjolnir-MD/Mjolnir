#ifndef MJOLNIR_STOICHIOMETRIC_PARAMETER_LIST_HPP
#define MJOLNIR_STOICHIOMETRIC_PARAMETER_LIST_HPP
#include <mjolnir/util/empty.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/forcefield/global/ParameterList.hpp>
#include <algorithm>

namespace mjolnir
{
// ----------------------------------------------------------------------------
// parameter combination rule for stoihchiometric interaction
// Now stoichiometric interaction only support uniform potential,
// so potentialT::parameter_type is empty_t.
//
template<typename traitsT, typename potentialT>
class StoichiometricEmptyCombinationRule final
    : public ParameterListBase<traitsT, potentialT>
{
  public:
    using base_type           = ParameterListBase<traitsT, potentialT>;
    using pair_parameter_type = typename potentialT::parameter_type;
    using real_type           = typename base_type::real_type;
    using parameter_type      = typename potentialT::parameter_type;

    // topology stuff
    using system_type          = System<traitsT>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traitsT>;

  public:
    StoichiometricEmptyCombinationRule(
        const std::vector<std::size_t>&& participants_a,
        const std::vector<std::size_t>&& participants_b,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : participants_(participants_a),
        participants_a_(participants_a), participants_b_(participants_b),
        exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        MJOLNIR_GET_DEFAULT_LOGGER_DEBUG();

        this->participants_.insert(participants_.end(),
                participants_b.begin(), participants_b.end());

        MJOLNIR_LOG_DEBUG("participants size is ", participants_.size());
        MJOLNIR_LOG_DEBUG("participants a size is ", participants_a_.size());
        MJOLNIR_LOG_DEBUG("participants b size is ", participants_b_.size());
    }

    StoichiometricEmptyCombinationRule(const StoichiometricEmptyCombinationRule&) = default;
    StoichiometricEmptyCombinationRule(StoichiometricEmptyCombinationRule&&) = default;
    StoichiometricEmptyCombinationRule& operator=(
            const StoichiometricEmptyCombinationRule&) = default;
    StoichiometricEmptyCombinationRule& operator=(
            StoichiometricEmptyCombinationRule&&) = default;
    ~StoichiometricEmptyCombinationRule() override {}

    pair_parameter_type prepare_params(std::size_t, std::size_t) const noexcept override
    {
        return empty_t{};
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potentialT& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol, pot);
        is_following_participants_.resize(sys.size(), false);
        for(const auto& i : participants_b_)
        {
            is_following_participants_[i] = true;
        }
        return;
    }

    void update(const system_type& sys, const topology_type& topol,
                const potentialT& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->max_cutoff_length_ = pot.max_cutoff();
        this->exclusion_list_.make(sys, topol);
        return;
    }

    real_type max_cutoff_length() const noexcept override {return max_cutoff_length_;}

    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        // i is leading_participant and j is a member of possible_partners_of
        // from implementation of SpatialPartition
        if(is_following_participants_[j] && i != j)
        {
            return !exclusion_list_.is_excluded(i, j);
        }
        return false;
    }
    std::vector<std::size_t> const& participants() const noexcept override
    {
        return this->participants_;
    }
    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return make_range(participants_a_.begin(), participants_a_.end());
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t, const std::size_t) const noexcept override
    {
        return make_range(participants_b_.begin(), participants_b_.end());
    }

    // used in GlobalStoichiometricInteraction
    std::vector<std::size_t> const& participants_a() const noexcept
    {
        return this->participants_a_;
    }
    std::vector<std::size_t> const& participants_b() const noexcept
    {
        return this->participants_b_;
    }

    base_type* clone() const override
    {
        return new StoichiometricEmptyCombinationRule(*this);
    }

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

  private:

    real_type                max_cutoff_length_;
    std::vector<std::size_t> participants_;
    std::vector<std::size_t> participants_a_;
    std::vector<std::size_t> participants_b_;
    std::vector<bool>        is_following_participants_;
    ExclusionList<traitsT>   exclusion_list_;
};

} // mjolnir
#endif // MJOLNIR_STOICHIOMETRIC_PARAMETER_LIST_HPP
