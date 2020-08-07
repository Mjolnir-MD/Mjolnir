#ifndef MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
#define MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
#include <mjolnir/forcefield/MultipleBasin/MultipleBasinUnitBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/ForceFieldBase.hpp>
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/core/GlobalForceField.hpp>
#include <mjolnir/core/ExternalForceField.hpp>
#include <mjolnir/util/string.hpp>
#include <algorithm>
#include <numeric>
#include <memory>

namespace mjolnir
{

// MultipleBasinForceField
//
// In some cases, protein has several domains that undergoes conformational
// changes independently.
//
//  .-.  .-.      __  .-.
// ( A )( B ) -> (A'|( B )
//  `-'  `-'      `-  `-'
//      |            |
//      v            v
//  .-.  __        __ __
// ( A )|B')  ->  (A'|B')
//  `-'  -'        `- -'
//
// To handle such cases, MultipleBasinForceField has several N-basin forcefields
// inside and manage them in a proper way.
//
template<typename traitsT>
class MultipleBasinForceField : public ForceFieldBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = ForceFieldBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = Topology;

    using local_forcefield_type      = LocalForceField<traits_type>;
    using global_forcefield_type     = GlobalForceField<traits_type>;
    using external_forcefield_type   = ExternalForceField<traits_type>;
    using constraint_forcefield_type = ConstraintForceField<traits_type>;
    using forcefield_type            = std::tuple<
        local_forcefield_type, global_forcefield_type, external_forcefield_type>;

    // a set of potentials that are correlated in the way of MultipleBasin.
    using multiple_basin_unit_type =
        std::unique_ptr<MultipleBasinUnitBase<traits_type>>;

  public:

    MultipleBasinForceField(forcefield_type&& common,
                            constraint_forcefield_type&& constraint,
                            std::vector<multiple_basin_unit_type>&& units)
        : loc_common_(std::move(std::get<0>(common))),
          glo_common_(std::move(std::get<1>(common))),
          ext_common_(std::move(std::get<2>(common))),
          constraint_(std::move(constraint)),
          units_(std::move(units))
    {}
    MultipleBasinForceField(forcefield_type&& common, // no constraint
                            std::vector<multiple_basin_unit_type>&& units)
        : loc_common_(std::move(std::get<0>(common))),
          glo_common_(std::move(std::get<1>(common))),
          ext_common_(std::move(std::get<2>(common))),
          units_(std::move(units))
    {}

    ~MultipleBasinForceField() override = default;
    MultipleBasinForceField(const MultipleBasinForceField&) = delete;
    MultipleBasinForceField(MultipleBasinForceField&&)      = default;
    MultipleBasinForceField& operator=(const MultipleBasinForceField&) = delete;
    MultipleBasinForceField& operator=(MultipleBasinForceField&&)      = default;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        // -------------------------------------------------------------------
        // write topologies. Here we assume topologies are the same

        MJOLNIR_LOG_INFO("writing topology");
        topol_.resize(sys.size());

        for(auto& unit : this->units_)
        {
            // each unit does not know the size of the system.
            // we need to tell them the size.
            unit->write_topology(sys, topol_);
        }
        loc_common_.write_topology(topol_);

        topol_.construct_molecules();

        MJOLNIR_LOG_INFO("initializing forcefields");
        for(auto& unit : this->units_)
        {
            unit->initialize(sys, this->topol_);
        }
        loc_common_.initialize(sys);
        glo_common_.initialize(sys, this->topol_);
        ext_common_.initialize(sys);
        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        for(const auto& unit : this->units_)
        {
            unit->calc_force(sys);
        }
        sys.preprocess_forces();
        loc_common_.calc_force(sys);
        glo_common_.calc_force(sys);
        ext_common_.calc_force(sys);
        sys.postprocess_forces();
        return ;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E = real_type(0.0);

        for(const auto& unit : this->units_)
        {
            E += unit->calc_energy(sys);
        }
        E += loc_common_.calc_energy(sys);
        E += glo_common_.calc_energy(sys);
        E += ext_common_.calc_energy(sys);
        return E;
    }

    void update(const system_type& sys) override
    {
        // update parameters (e.g. temperature). TODO: topologies?

        for(auto& unit : this->units_)
        {
            unit->update(sys, this->topol_);
        }
        loc_common_.update(sys);
        glo_common_.update(sys, this->topol_);
        ext_common_.update(sys);
        return;
    }

    // update margin of neighbor list
    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        for(auto& unit : this->units_)
        {
            unit->reduce_margin(dmargin, sys);
        }
        loc_common_.reduce_margin(dmargin, sys);
        glo_common_.reduce_margin(dmargin, sys);
        ext_common_.reduce_margin(dmargin, sys);
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        for(auto& unit : this->units_)
        {
            unit->scale_margin(scale, sys);
        }
        loc_common_.scale_margin(scale, sys);
        glo_common_.scale_margin(scale, sys);
        ext_common_.scale_margin(scale, sys);
        return;
    }

    // -----------------------------------------------------------------------
    // energy output format

    void format_energy_name(std::string& fmt) const override
    {
        using namespace mjolnir::literals::string_literals;

        for(const auto& unit : this->units_)
        {
            unit->format_energy_name(fmt);
        }

        fmt += "Common{"_s;
        loc_common_.format_energy_name(fmt);
        glo_common_.format_energy_name(fmt);
        ext_common_.format_energy_name(fmt);
        fmt += "} "_s;
        return;
    }

    real_type format_energy(const system_type& sys, std::string& fmt) const override
    {
        real_type total = 0.0;
        for(const auto& unit : this->units_)
        {
            total += unit->format_energy(sys, fmt);
        }

        fmt += "Common{"_s;
        total += loc_common_.format_energy(sys, fmt);
        total += glo_common_.format_energy(sys, fmt);
        total += ext_common_.format_energy(sys, fmt);
        fmt += "} "_s;
        return total;
    }

    // -----------------------------------------------------------------------

    topology_type const&              topology()   const noexcept override {return topol_;}
    constraint_forcefield_type const& constraint() const noexcept override {return constraint_;}

    local_forcefield_type    const& common_local()    const noexcept {return loc_common_;}
    global_forcefield_type   const& common_global()   const noexcept {return glo_common_;}
    external_forcefield_type const& common_external() const noexcept {return ext_common_;}

    std::vector<multiple_basin_unit_type> const& units() const noexcept {return units_;}

  private:

    topology_type            topol_;
    local_forcefield_type    loc_common_;
    global_forcefield_type   glo_common_;
    external_forcefield_type ext_common_;

    constraint_forcefield_type constraint_;

    std::vector<multiple_basin_unit_type> units_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class MultipleBasinForceField<SimulatorTraits<double, UnlimitedBoundary       >>;
extern template class MultipleBasinForceField<SimulatorTraits<float,  UnlimitedBoundary       >>;
extern template class MultipleBasinForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class MultipleBasinForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
