#ifndef MJOLNIR_TEST_UTIL_STUB_POTENTIAL_HPP
#define MJOLNIR_TEST_UTIL_STUB_POTENTIAL_HPP
#include <mjolnir/math/Vector.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/forcefield/global/ParameterList.hpp>

namespace mjolnir
{
namespace test
{

template<typename realT>
class StubPotential
{
  public:
    using real_type = realT;
    struct parameter_type
    {
        std::size_t i;
        std::size_t j;
    };

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.0);
    }

  public:

    StubPotential() noexcept {}
    ~StubPotential() = default;

    real_type potential(const real_type, const parameter_type&) const noexcept
    {
        return 0.0;
    }
    real_type derivative(const real_type, const parameter_type&) const noexcept
    {
        return 0.0;
    }

    template<typename T>
    void initialize(const System<T>&) noexcept {return;}

    template<typename T>
    void update(const System<T>&) noexcept {return;}

    template<typename InputIterator>
    real_type max_cutoff(const InputIterator, const InputIterator) const noexcept
    {
        return 2.0;
    }
    real_type absolute_cutoff(const parameter_type&) const noexcept
    {
        return 2.0;
    }

    static const char* name() noexcept {return "Stub";}
};

template<typename traitsT>
struct StubParameterList
    : public ParameterListBase<traitsT, StubPotential<typename traitsT::real_type>>
{
    using traits_type          = traitsT;
    using real_type            = typename traits_type::real_type;
    using potential_type       = StubPotential<real_type>;
    using base_type            = ParameterListBase<traits_type, potential_type>;

    using pair_parameter_type  = typename potential_type::parameter_type;

    using system_type          = System<traits_type>;
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using group_id_type        = typename topology_type::group_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;
    using ignore_group_type    = IgnoreGroup   <group_id_type>;
    using exclusion_list_type  = ExclusionList<traits_type>;

    explicit StubParameterList(const double cutoff,
        const std::vector<std::size_t>& parameters,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_mol, ignore_group_type ignore_grp)
        : cutoff_(cutoff), participants_(parameters),
          exclusion_list_(exclusions, std::move(ignore_mol), std::move(ignore_grp))
    {}

    real_type max_cutoff_length() const noexcept override {return cutoff_;}

    pair_parameter_type prepare_params(std::size_t i, std::size_t j) const noexcept override
    {
        return pair_parameter_type{i, j};
    }

    void initialize(const system_type& sys, const topology_type& topol,
                    const potential_type& pot) noexcept override
    {
        this->update(sys, topol, pot);
        return;
    }
    void update(const system_type& sys, const topology_type& topol,
                const potential_type&) noexcept override
    {
        // update exclusion list based on sys.topology()
        exclusion_list_.make(sys, topol);
        return;
    }

    std::vector<std::size_t> const& participants() const noexcept override
    {
        return this->participants_;
    }
    mjolnir::range<typename std::vector<std::size_t>::const_iterator>
    leading_participants() const noexcept override
    {
        return mjolnir::make_range(participants_.begin(), std::prev(participants_.end()));
    }
    mjolnir::range<typename std::vector<std::size_t>::const_iterator>
    possible_partners_of(const std::size_t participant_idx,
                         const std::size_t /*particle_idx*/) const noexcept override
    {
        return mjolnir::make_range(participants_.begin() + participant_idx + 1,
                                   participants_.end());
    }
    bool has_interaction(const std::size_t i, const std::size_t j) const noexcept override
    {
        return (i < j);
    }

    std::string name() const {return "Stub";}

    exclusion_list_type const& exclusion_list() const noexcept override
    {
        return exclusion_list_; // for testing
    }

    base_type* clone() const override
    {
        return new StubParameterList(*this);
    }

  private:

    real_type cutoff_;
    std::vector<std::size_t> participants_;
    exclusion_list_type  exclusion_list_;
};

} // test
} // mjolnir
#endif// MJOLNIR_TEST_UTIL_CHECK_FORCE_HPP
