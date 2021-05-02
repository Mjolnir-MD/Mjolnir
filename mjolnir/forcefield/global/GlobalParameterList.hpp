#ifndef MJOLNIR_GLOBAL_PARAMETER_LIST_HPP
#define MJOLNIR_GLOBAL_PARAMETER_LIST_HPP
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <vector>
#include <map>

namespace mjolnir
{

template<typename traitsT, typename potentialT>
struct CombinationRuleBase
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using pair_parameter_type = typename potential_type::parameter_type;

  public:

    virtual ~CombinationRuleBase() = default;

    // calculate combination of parameters
    virtual pair_parameter_type
    generate(const std::size_t i, const std::size_t j) const noexcept = 0;

    // This function is for cutoff distance
    virtual real_type max_cutoff_length(const potential_type& pot) const = 0;

    virtual std::vector<std::size_t>&       participants()       noexcept = 0;
    virtual std::vector<std::size_t> const& participants() const noexcept = 0;
};

//
// GlobalParameterList contains list of per-particle parameter used in a pair
// interaction and a combination rule.
//
template<typename traitsT, typename potentialT>
class GlobalParameterList
{
  public:
    using traits_type           = traitsT;
    using potential_type        = potentialT;
    using combination_rule_base = CombinationRuleBase<traits_type, potential_type>;
    using combination_rule_type = std::unique_ptr<combination_rule_base>;
    using pair_parameter_type   = typename potential_type::parameter_type;
    using system_type           = System<traits_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

  public:

    GlobalParameterList(combination_rule_type rule,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : rule_(std::move(rule)),
        exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {}

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept
    {
        return rule_->generate(i, j);
    }

    // These value

    template<typename PotentialT>
    void initialize(const system_type& sys, const topology_type& topol,
                    const PotentialT& pot) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol, pot);
        return;
    }

    template<typename PotentialT>
    void update(const system_type& sys, const topology_type& topol,
                const PotentialT& pot) noexcept
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->max_cutoff_length_ = rule_->max_cutoff_length(pot);
        this->exclusion_list_.make(sys, topol);
        return;
    }

    real_type max_cutoff_length() const noexcept {return max_cutoff_length_;}

    // -----------------------------------------------------------------------
    // Here, the default implementation uses Newton's 3rd law to reduce
    // calculation. For an interacting pair (i, j), forces applied to i and j
    // are equal in magnitude and opposite in direction. So, if a pair (i, j) is
    // listed, (j, i) is not needed.
    //     See implementation of VerletList, CellList and GlobalPairInteraction
    // for more details about the usage of these functions.
    //
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    std::vector<std::size_t> const& participants() const noexcept
    {
        return rule_->participants();
    }
    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept
    {
        const auto& ps = rule_->participants();
        return make_range(ps.begin(), std::prev(ps.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept
    {
        return make_range(ps.begin() + participant_idx + 1, ps.end());
    }

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    exclusion_list_type const& exclusion_list() const noexcept
    {
        return exclusion_list_;
    }

    combination_rule_type&       rule()       noexcept {return rule_;}
    combination_rule_type const& rule() const noexcept {return rule_;}

  private:

    real_type                max_cutoff_length_;
    combination_rule_type    rule_;
    exclusion_list_type      exclusion_list_;
};

// ----------------------------------------------------------------------------
// widely-used parameter combination rules
//
// CombinationRules should provide the following two functionality.
// 1. combination rule
// 2. the pair parameter that provides the maximum cutoff distance
//
// The second one depends on the functional form of the potential, so the
// comparator used to find the max cutoff distance is given by the potential
// class via get_parameter_comparator().
//     To avoid O(n^2) search while finding the maximum cutoff distance, it
// linearly search the parameters using the comparator. It means that, the
// combination of the maximum parameters (two identical, maximum parameters)
// should give the maximum cutoff distance. This is automatically satisfied
// if you are using Lorentz-Berthelot rule or some other variant of it. If
// you are using combination table, the comparator directly compare the
// pair parameters. So, the comparator should be able to be used for both
// per-particle parameters and pair parameters.

template<typename traitsT, typename potentialT>
class LorentzBerthelotRule final
    : public CombinationRuleBase<traitsT, potentialT>
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using real_type           = typename traits_type::real_type;
    using parameter_type      = std::pair<real_type, real_type>;
    using pair_parameter_type = std::pair<real_type, real_type>;
    using container_type      = std::vector<parameter_type>;

  public:

    ~LorentzBerthelotRule() override {}

    explicit LorentzBerthelotRule(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters)
    {
        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(this->parameters_.size() <= idx)
            {
                this->parameters_.resize(idx+1, default_parameter);
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }

    // calculate combination of parameters
    pair_parameter_type
    generate(const std::size_t i, const std::size_t j) const noexcept override
    {
        const auto& para1 = parameters_[i];
        const auto& para2 = parameters_[j];

        const auto sgm1 = para1.first;
        const auto eps1 = para1.second;
        const auto sgm2 = para2.first;
        const auto eps2 = para2.second;

        return std::make_pair((sgm1 + sgm2) / 2,
                             ((eps1 == eps2) ? eps1 : std::sqrt(eps1 * eps2)));
    }

    // This function is for cutoff distance
    real_type max_cutoff_length(const potential_type& pot) const override
    {
        if(parameters_.empty())
        {
            return std::numeric_limits<real_type>::infinity();
        }
        const auto max_iter = std::max_element(
                parameters_.begin(), parameters_.end(), pot.get_parameter_comparator());

        return PotentialT(*max_iter).cutoff();
    }

    // accessors for testing

    container_type&       parameters()       noexcept {return parameters_;}
    container_type const& parameters() const noexcept {return parameters_;}

    std::vector<std::size_t>&       participants()       noexcept override {return participants_;}
    std::vector<std::size_t> const& participants() const noexcept override {return participants_;}

  private:

    std::vector<std::size_t> participants_;
    container_type           parameters_;
};

template<typename traitsT, typename potentialT>
struct CombinationTable
    : public CombinationRuleBase<traitsT, potentialT>
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using parameter_type      = std::string;
    using pair_parameter_type = typename potential_type::parameter_type;
    using container_type      = std::vector<parameter_type>;
    using table_type          = std::unordered_map<std::string, pair_parameter_type>;

  public:

    ~CombinationTable() override {}

    CombinationTable(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        table_type&& table)
        : table_(std::move(table))
    {
        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(this->parameters_.size() <= idx)
            {
                this->parameters_.resize(idx+1, default_parameter);
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }

    // Since it concats the names, it is better to keep parameter name short.
    pair_parameter_type
    generate(const std::size_t i, const std::size_t j) const noexcept override
    {
        const auto& para1 = parameters_[i];
        const auto& para2 = parameters_[j];
        return table_[para1 + para2];
    }

    // This function is for cutoff distance
    real_type max_cutoff_length(const potential_type& pot) const override
    {
        if(parameters_.empty())
        {
            return std::numeric_limits<real_type>::infinity();
        }

        using value_type = typename table_type::value_type; // kv-pair
        const auto max_iter = std::max_element(table_.begin(), table_.end(),
            [&](const value_type& lkv, const value_type& rkv) noexcept -> bool {
                return pot.get_parameter_comparator(lkv.second, rkv.second);
            });

        return PotentialT(max_iter->second).cutoff();
    }

    // accessors

    table_type const& table() const noexcept {return table_;}
    table_type&       table()       noexcept {return table_;}

    container_type&       parameters()       noexcept {return parameters_;}
    container_type const& parameters() const noexcept {return parameters_;}

    std::vector<std::size_t>&       participants()       noexcept override {return participants_;}
    std::vector<std::size_t> const& participants() const noexcept override {return participants_;}

  private:

    table_type               table_;
    std::vector<std::size_t> participants_;
    container_type           parameters_;
};

} // mjolnir
#endif// MJOLNIR_GLOBAL_PARAMETER_LIST_HPP
