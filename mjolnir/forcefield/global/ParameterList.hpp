#ifndef MJOLNIR_GLOBAL_PARAMETER_LIST_HPP
#define MJOLNIR_GLOBAL_PARAMETER_LIST_HPP
#include <mjolnir/util/empty.hpp>
#include <mjolnir/core/ExclusionList.hpp>
#include <mjolnir/core/System.hpp>
#include <vector>
#include <map>

namespace mjolnir
{
//
// ParameterList is responsible for generating pair-parameter for global
// interaction. It can be decomposed into the following roles.
//
// - Store per-particle parameters
// - Generate pair-parameters using combination rule (like Lorentz-Berthelot)
// - Store topology-dependent interaction rule (exclusion list)
//
// It requires potentialT to have
// - typedef potentialT::parameter_type
// -         potentialT::potentialT(parameter_type) (used in SpatialPartition)
//
// The default implementation requires
// - typedef    potentialT::parameter_comparator
// - non-static potentialT::cutoff()
//
template<typename traitsT, typename potentialT>
struct ParameterListBase
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using real_type           = typename traits_type::real_type;
    using pair_parameter_type = typename potential_type::parameter_type;
    using system_type         = System<traits_type>;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

  public:

    virtual ~ParameterListBase() = default;

    virtual void initialize(const system_type&, const topology_type&, const potential_type&) noexcept = 0;
    virtual void update    (const system_type&, const topology_type&, const potential_type&) noexcept = 0;

    // calculate combination of parameters
    virtual pair_parameter_type
    prepare_params(const std::size_t i, const std::size_t j) const noexcept = 0;

    virtual real_type max_cutoff_length() const noexcept = 0;

    // Topology-related filtering
    virtual bool has_interaction(const std::size_t i, const std::size_t j) const noexcept = 0;

    virtual std::vector<std::size_t> const& participants() const noexcept = 0;

    virtual range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept = 0;

    virtual range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t particle_idx) const noexcept = 0;

    virtual exclusion_list_type const& exclusion_list() const noexcept = 0;

    virtual ParameterListBase* clone() const = 0;
};

// This class just contains unique_ptr<ParameterListBase> and forwards functions
// into Base*. This makes parameter list handling a bit easier.
template<typename traitsT, typename potentialT>
struct ParameterList
{
  public:
    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using base_type           = ParameterListBase<traits_type, potential_type>;
    using real_type           = typename base_type::real_type;
    using pair_parameter_type = typename base_type::pair_parameter_type;

    // topology stuff
    using system_type          = System<traits_type>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

  public:

    explicit ParameterList(std::unique_ptr<base_type> ptr)
       : ptr_(std::move(ptr))
    {}
    ~ParameterList() = default;

    ParameterList(ParameterList const& other)
        : ptr_(other.ptr_->clone())
    {}
    ParameterList& operator=(ParameterList const& other)
    {
        this->ptr_.reset(other.ptr_->clone());
        return *this;
    }
    ParameterList(ParameterList&&)            = default;
    ParameterList& operator=(ParameterList&&) = default;

    base_type&        ref() const noexcept {return *ptr_;}
    base_type const& cref() const noexcept {return *ptr_;}

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept
    {
        return ptr_->initialize(sys, topol, pot);
    }
    void update    (const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept
    {
        return ptr_->initialize(sys, topol, pot);
    }

    // calculate combination of parameters
    pair_parameter_type
    prepare_params(const std::size_t i, const std::size_t j) const noexcept
    {
        return ptr_->prepare_params(i, j);
    }

    real_type max_cutoff_length() const noexcept
    {
        return ptr_->max_cutoff_length();
    }

    // Topology-related filtering
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept
    {
        return ptr_->has_interaction(i, j);
    }

    std::vector<std::size_t> const& participants() const noexcept
    {
        return ptr_->participants();
    }

    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept
    {
        return ptr_->leading_participants();
    }

    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t particle_idx) const noexcept
    {
        return ptr_->possible_partners_of(participant_idx, particle_idx);
    }

    exclusion_list_type const& exclusion_list() const noexcept
    {
        return ptr_->exclusion_list();
    }

  private:

    std::unique_ptr<ParameterListBase<traits_type, potential_type>> ptr_;
};

// ----------------------------------------------------------------------------
// widely-used parameter combination rules

// Lorentz-Berthelot Rule uses sigma and epsilon and
// pair-parameters are defined as
//
// sigma_ij   = (sigma_i + sigma_j) / 2
// epsilon_ij = sqrt(epsilon_i * epsilon_j)
//
template<typename traitsT, typename potentialT>
class LorentzBerthelotRule final
    : public ParameterListBase<traitsT, potentialT>
{
  public:
    using traits_type           = traitsT;
    using potential_type        = potentialT;
    using base_type             = ParameterListBase<traits_type, potential_type>;

    // Lorentz-Bertherot rule requires the parameter to have sigma and epsilon.
    // The `parameter_type`s should have `sigma` and `epsilon` as a member var
    // of struct.
    using real_type             = typename traits_type::real_type;
    using parameter_type        = typename potential_type::parameter_type;
    using pair_parameter_type   = typename potential_type::parameter_type;
    using container_type        = std::vector<parameter_type>;

    // topology stuff
    using system_type           = System<traits_type>;
    using topology_type         = Topology;
    using molecule_id_type      = typename topology_type::molecule_id_type;
    using group_id_type         = typename topology_type::group_id_type;
    using connection_kind_type  = typename topology_type::connection_kind_type;
    using ignore_molecule_type  = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type     = IgnoreGroup   <group_id_type>;
    using exclusion_list_type   = ExclusionList <traits_type>;

  public:

    LorentzBerthelotRule(
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        assert(parameters.size() != 0);

        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(this->parameters_.size() <= idx)
            {
                this->parameters_.resize(idx+1, {0, 0});
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }

    LorentzBerthelotRule(const LorentzBerthelotRule&) = default;
    LorentzBerthelotRule(LorentzBerthelotRule&&) = default;
    LorentzBerthelotRule& operator=(const LorentzBerthelotRule&) = default;
    LorentzBerthelotRule& operator=(LorentzBerthelotRule&&) = default;
    ~LorentzBerthelotRule() override {}

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        const auto& para1 = parameters_[i];
        const auto& para2 = parameters_[j];

        const auto sgm1 = para1.sigma;
        const auto eps1 = para1.epsilon;
        const auto sgm2 = para2.sigma;
        const auto eps2 = para2.epsilon;

        return pair_parameter_type{(sgm1 + sgm2) / 2,
                                   (eps1 == eps2) ? eps1 : std::sqrt(eps1 * eps2)};
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol, pot);
        return;
    }

    void update(const system_type& sys, const topology_type& topol,
                const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->max_cutoff_length_ = pot.max_cutoff(parameters_.begin(), parameters_.end());
        this->exclusion_list_.make(sys, topol);
        return;
    }

    real_type max_cutoff_length() const noexcept override {return max_cutoff_length_;}

    // -----------------------------------------------------------------------
    // Here, the default implementation uses Newton's 3rd law to reduce
    // calculation. For an interacting pair (i, j), forces applied to i and j
    // are equal in magnitude and opposite in direction. So, if a pair (i, j) is
    // listed, (j, i) is not needed.
    //     See implementation of VerletList, CellList and GlobalPairInteraction
    // for more details about the usage of these functions.
    //
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    std::vector<std::size_t> const& participants() const noexcept override
    {
        return this->participants_;
    }
    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept override
    {
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }

    base_type* clone() const override
    {
        return new LorentzBerthelotRule(*this);
    }

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    container_type&       parameters()       noexcept {return parameters_;}
    container_type const& parameters() const noexcept {return parameters_;}

  private:

    real_type                max_cutoff_length_;
    std::vector<std::size_t> participants_;
    container_type           parameters_;
    exclusion_list_type      exclusion_list_;
};

// =============================================================================
//
// Combination Table contains a map from "name:name" to parameter and requires
// each particle to have its name.
//
template<typename traitsT, typename potentialT>
struct CombinationTable final
    : public ParameterListBase<traitsT, potentialT>
{
  public:

    using traits_type         = traitsT;
    using potential_type      = potentialT;
    using base_type           = ParameterListBase<traits_type, potential_type>;

    using real_type           = typename traits_type::real_type;
    using parameter_type      = std::string;
    using pair_parameter_type = typename potential_type::parameter_type;
    using container_type      = std::vector<parameter_type>;
    using table_type          = std::unordered_map<std::string, pair_parameter_type>;

    // topology stuff
    using system_type          = System<traits_type>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList <traits_type>;

  public:

    CombinationTable(table_type&& table,
        const std::vector<std::pair<std::size_t, parameter_type>>& parameters,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : table_(std::move(table)),
        exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        assert(parameters.size() != 0);

        this->parameters_  .reserve(parameters.size());
        this->participants_.reserve(parameters.size());
        for(const auto& idxp : parameters)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(this->parameters_.size() <= idx)
            {
                this->parameters_.resize(idx+1, std::string(""));
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    CombinationTable(const CombinationTable&) = default;
    CombinationTable(CombinationTable&&)      = default;
    CombinationTable& operator=(const CombinationTable&) = default;
    CombinationTable& operator=(CombinationTable&&)      = default;
    ~CombinationTable() override {}

    // Since it concats the names, it is better to keep parameter name short.
    pair_parameter_type
    prepare_params(const std::size_t i, const std::size_t j) const noexcept override
    {
        const auto& para1 = parameters_[i];
        const auto& para2 = parameters_[j];
        return table_.at(para1 + std::string(":") + para2);
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->update(sys, topol, pot);
        return;
    }

    void update(const system_type& sys, const topology_type& topol,
                const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(parameters_.empty()) // set a stub
        {
            this->max_cutoff_length_ = 1.0;
        }
        else // find max_cutoff from pair-parameters in the table
        {
            this->max_cutoff_length_ = 0.0;
            for(const auto& kv : this->table_)
            {
                const auto& parameter = kv.second;

                this->max_cutoff_length_ = std::max(this->max_cutoff_length_,
                        pot.absolute_cutoff(parameter));
            }
        }
        this->exclusion_list_.make(sys, topol);
        return;
    }

    real_type max_cutoff_length() const noexcept override {return max_cutoff_length_;}

    // -----------------------------------------------------------------------
    // Here, the default implementation uses Newton's 3rd law to reduce
    // calculation. For an interacting pair (i, j), forces applied to i and j
    // are equal in magnitude and opposite in direction. So, if a pair (i, j) is
    // listed, (j, i) is not needed.
    //     See implementation of VerletList, CellList and GlobalPairInteraction
    // for more details about the usage of these functions.
    //
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    std::vector<std::size_t> const& participants() const noexcept override
    {
        return this->participants_;
    }
    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept override
    {
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }

    base_type* clone() const override
    {
        return new CombinationTable(*this);
    }

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    table_type const& table() const noexcept {return table_;}
    table_type&       table()       noexcept {return table_;}

    container_type&       parameters()       noexcept {return parameters_;}
    container_type const& parameters() const noexcept {return parameters_;}

  private:

    real_type                max_cutoff_length_;
    table_type               table_;
    std::vector<std::size_t> participants_;
    container_type           parameters_;
    exclusion_list_type      exclusion_list_;
};

// ----------------------------------------------------------------------------
// for potentials that do not have any parameter.
// The potential::parameter_type should be mjolnir::empty_t.

template<typename traitsT, typename potentialT>
class EmptyCombinationRule final
    : public ParameterListBase<traitsT, potentialT>
{
  public:
    using traits_type           = traitsT;
    using potential_type        = potentialT;
    using base_type             = ParameterListBase<traits_type, potential_type>;
    using real_type             = typename traits_type::real_type;
    using pair_parameter_type   = typename potential_type::parameter_type;

    static_assert(std::is_same<pair_parameter_type, empty_t>::value, "");

    // topology stuff
    using system_type           = System<traits_type>;
    using topology_type         = Topology;
    using molecule_id_type      = typename topology_type::molecule_id_type;
    using group_id_type         = typename topology_type::group_id_type;
    using connection_kind_type  = typename topology_type::connection_kind_type;
    using ignore_molecule_type  = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type     = IgnoreGroup   <group_id_type>;
    using exclusion_list_type   = ExclusionList <traits_type>;

  public:

    EmptyCombinationRule(
        const std::vector<std::size_t>& participants,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : participants_should_be_updated_(false), participants_(participants),
        exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {
        assert(participants.size() != 0);
    }
    EmptyCombinationRule(
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
      : participants_should_be_updated_ (true),
        exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {}

    EmptyCombinationRule(const EmptyCombinationRule&) = default;
    EmptyCombinationRule(EmptyCombinationRule&&) = default;
    EmptyCombinationRule& operator=(const EmptyCombinationRule&) = default;
    EmptyCombinationRule& operator=(EmptyCombinationRule&&) = default;
    ~EmptyCombinationRule() override {}

    pair_parameter_type prepare_params(std::size_t, std::size_t) const noexcept override
    {
        return empty_t{};
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        if(participants_should_be_updated_)
        {
            participants_.resize(sys.size());
            std::iota(participants_.begin(), participants_.end(), 0);

            participants_should_be_updated_ = false;
        }
        this->update(sys, topol, pot);
        return;
    }

    void update(const system_type& sys, const topology_type& topol,
                const potential_type& pot) noexcept override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->max_cutoff_length_ = pot.absolute_cutoff(empty_t{});
        this->exclusion_list_.make(sys, topol);
        return;
    }

    real_type max_cutoff_length() const noexcept override {return max_cutoff_length_;}

    // -----------------------------------------------------------------------
    // Here, the default implementation uses Newton's 3rd law to reduce
    // calculation. For an interacting pair (i, j), forces applied to i and j
    // are equal in magnitude and opposite in direction. So, if a pair (i, j) is
    // listed, (j, i) is not needed.
    //     See implementation of VerletList, CellList and GlobalPairInteraction
    // for more details about the usage of these functions.
    //
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        // if not excluded, the pair has interaction.
        return (i < j) && !exclusion_list_.is_excluded(i, j);
    }
    std::vector<std::size_t> const& participants() const noexcept override
    {
        return this->participants_;
    }
    range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return make_range(participants_.begin(), std::prev(participants_.end()));
    }
    range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept override
    {
        return make_range(participants_.begin() + participant_idx + 1,
                          participants_.end());
    }

    base_type* clone() const override
    {
        return new EmptyCombinationRule(*this);
    }

    // ------------------------------------------------------------------------
    // the following accessers would be used in tests.

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_;
    }

  private:

    bool                     participants_should_be_updated_;
    real_type                max_cutoff_length_;
    std::vector<std::size_t> participants_;
    exclusion_list_type      exclusion_list_;
};

} // mjolnir
#endif// MJOLNIR_GLOBAL_PARAMETER_LIST_HPP
