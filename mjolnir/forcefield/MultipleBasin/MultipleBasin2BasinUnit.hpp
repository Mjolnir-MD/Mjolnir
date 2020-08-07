#ifndef MJOLNIR_CORE_MULTIPLE_BASIN_2_BASIN_FORCE_FIELD_HPP
#define MJOLNIR_CORE_MULTIPLE_BASIN_2_BASIN_FORCE_FIELD_HPP
#include <mjolnir/forcefield/MultipleBasin/MultipleBasinUnitBase.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/LocalForceField.hpp>
#include <mjolnir/core/GlobalForceField.hpp>
#include <mjolnir/core/ExternalForceField.hpp>
#include <mjolnir/util/string.hpp>
#include <algorithm>
#include <numeric>
#include <memory>

namespace mjolnir
{

// 2-basin MultipleBasin unit.
//
template<typename traitsT>
class MultipleBasin2BasinUnit final: public MultipleBasinUnitBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = MultipleBasinUnitBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using topology_type   = typename base_type::topology_type;
    using local_forcefield_type    = LocalForceField<traits_type>;
    using global_forcefield_type   = GlobalForceField<traits_type>;
    using external_forcefield_type = ExternalForceField<traits_type>;
    using forcefield_type          = std::tuple<
        local_forcefield_type, global_forcefield_type, external_forcefield_type>;
    using coordinate_container_type =
        typename system_type::coordinate_container_type;

  public:

    MultipleBasin2BasinUnit(
            const real_type    delta,
            const std::string& name1, const std::string& name2,
            const real_type    dV1,   const real_type    dV2,
            forcefield_type&&  ff1,   forcefield_type&&  ff2)
        : dV1_(dV1), dV2_(dV2), delta_(delta),
          rdelta_(real_type(1.0) / delta), delta_sq_(delta * delta),
          c2_over_c1_(0.0), name1_(name1), name2_(name2),
          loc1_(std::move(std::get<0>(ff1))), loc2_(std::move(std::get<0>(ff2))),
          glo1_(std::move(std::get<1>(ff1))), glo2_(std::move(std::get<1>(ff2))),
          ext1_(std::move(std::get<2>(ff1))), ext2_(std::move(std::get<2>(ff2)))
    {}

    ~MultipleBasin2BasinUnit() override = default;
    MultipleBasin2BasinUnit(const MultipleBasin2BasinUnit&) = default;
    MultipleBasin2BasinUnit(MultipleBasin2BasinUnit&&)      = default;
    MultipleBasin2BasinUnit& operator=(const MultipleBasin2BasinUnit&) = default;
    MultipleBasin2BasinUnit& operator=(MultipleBasin2BasinUnit&&)      = default;

    void write_topology(const system_type& sys, topology_type& topol) const override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        MJOLNIR_LOG_INFO("checking topologies are the same");

        Topology topol1(sys.size()), topol2(sys.size());

        loc1_.write_topology(topol1);
        loc2_.write_topology(topol2);
        topol1.construct_molecules();
        topol2.construct_molecules();

        if(topol1 != topol2)
        {
            MJOLNIR_LOG_ERROR("topologies of 2 basins (", name1_, " and ",
                              name2_, ") are different from each other.");
            MJOLNIR_LOG_ERROR("MultipleBasin does not support such a case.");
            throw std::runtime_error("mjolnir::MultipleBasin2BasinUnit: "
                    "Topologies of 2 basins shouold be the same");
        }

        MJOLNIR_LOG_INFO("writing topology");
        loc1_.write_topology(topol);
        return;
    }

    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->force_buffer_.resize(sys.size(),
                math::make_coordinate<coordinate_type>(0, 0, 0));

        // -------------------------------------------------------------------
        MJOLNIR_LOG_INFO("initializing topology of the first basin");

        loc1_.initialize(sys);
        glo1_.initialize(sys, topol);
        ext1_.initialize(sys);

        // -------------------------------------------------------------------
        MJOLNIR_LOG_INFO("initializing topology of the first basin");

        loc2_.initialize(sys);
        glo2_.initialize(sys, topol);
        ext2_.initialize(sys);

        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        // -------------------------------------------------------------------
        // calc force of V_MB first.

        sys.preprocess_forces();
        const auto V_1 = this->calc_force_and_energy_basin1(sys) + dV1_;
        sys.postprocess_forces();
        {
            // save the current forces to force_buffer_.
            // force_buffer is zero-cleared (at the end of this function),
            // so the forces in the system will be zero-cleared after this.
            using std::swap;
            swap(this->force_buffer_, sys.forces());
        }
        sys.preprocess_forces();
        const auto V_2 = this->calc_force_and_energy_basin2(sys) + dV2_;
        sys.postprocess_forces();

        const auto V_diff = V_1 - V_2;
        const auto coef   = V_diff / std::sqrt(V_diff * V_diff + 4 * delta_sq_);

        // F_V_MB = 1/2 [1 - Vdiff / sqrt(Vdiff^2 + 4 delta^2)] * F_V_1 +
        //          1/2 [1 + Vdiff / sqrt(Vdiff^2 + 4 delta^2)] * F_V_2

        const auto coef1 = real_type(0.5) * (real_type(1) - coef);
        const auto coef2 = real_type(0.5) * (real_type(1) + coef);

        // here, sys.forces has forces of basin2. force_buffer has forces of V1.
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) *= coef2;
            sys.force(i) += coef1 * force_buffer_[i];
            force_buffer_[i] = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        return ;
    }
    real_type calc_energy(const system_type& sys) const noexcept override
    {
        // -------------------------------------------------------------------
        // calc energy of V_MB
        const auto V_1    = this->calc_energy_basin1(sys) + this->dV1_;
        const auto V_2    = this->calc_energy_basin2(sys) + this->dV2_;
        const auto V_diff = V_1 - V_2;

        const auto V_MB =
            (V_1 + V_2 - std::sqrt(V_diff * V_diff + 4 * delta_sq_)) / 2;

        // calc eigen vector
        this->c2_over_c1_ = (V_MB - V_1) * this->rdelta_;

        return V_MB;
    }

    void update(const system_type& sys, const topology_type& topol) override
    {
        // update parameters (e.g. temperature). TODO: topologies?

        loc1_.update(sys);
        glo1_.update(sys, topol);
        ext1_.update(sys);

        loc2_.update(sys);
        glo2_.update(sys, topol);
        ext2_.update(sys);
        return;
    }

    // update margin of neighbor list
    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        loc1_.reduce_margin(dmargin, sys);
        glo1_.reduce_margin(dmargin, sys);
        ext1_.reduce_margin(dmargin, sys);

        loc2_.reduce_margin(dmargin, sys);
        glo2_.reduce_margin(dmargin, sys);
        ext2_.reduce_margin(dmargin, sys);
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        loc1_.scale_margin(scale, sys);
        glo1_.scale_margin(scale, sys);
        ext1_.scale_margin(scale, sys);

        loc2_.scale_margin(scale, sys);
        glo2_.scale_margin(scale, sys);
        ext2_.scale_margin(scale, sys);
        return;
    }

    void format_energy_name(std::string& fmt) const override
    {
        using namespace mjolnir::literals::string_literals;
        // Basin1{BondLength:Harmonic ...} Basin2{...} V_MB chi

        fmt += name1_ + "{"_s;
        loc1_.format_energy_name(fmt);
        glo1_.format_energy_name(fmt);
        ext1_.format_energy_name(fmt);
        if(!fmt.empty() && fmt.back() == ' ') {fmt.pop_back();}
        fmt += "} "_s;

        fmt += name2_ + "{"_s;
        loc2_.format_energy_name(fmt);
        glo2_.format_energy_name(fmt);
        ext2_.format_energy_name(fmt);
        if(!fmt.empty() && fmt.back() == ' ') {fmt.pop_back();}
        fmt += "} "_s;

        fmt += "MultipleBasinTotal ";
        fmt += "MultipleBasinChi ";
        return;
    }

    real_type format_energy(const system_type& sys, std::string& fmt) const override
    {
        using namespace mjolnir::literals::string_literals;
        // Basin1{BondLength:Harmonic ...} Basin2{...} V_MB chi

        fmt += name1_ + "{"_s;
        real_type V_1 = this->dV1_;
        V_1 += loc1_.format_energy(sys, fmt);
        V_1 += glo1_.format_energy(sys, fmt);
        V_1 += ext1_.format_energy(sys, fmt);
        fmt += "} "_s;

        fmt += name2_ + "{"_s;
        real_type V_2 = this->dV2_;
        V_2 += loc2_.format_energy(sys, fmt);
        V_2 += glo2_.format_energy(sys, fmt);
        V_2 += ext2_.format_energy(sys, fmt);
        fmt += "} "_s;

        // ( V1-V_MB   delta    ) (c1) = (0)
        // (  delta  V2+dV-V_MB ) (c2)   (0)
        // <=> (V1-V_MB) c1 + delta c2 = 0
        // <=> (V_MB-V1) c1 = delta c2

        const auto V_diff =  V_1 - V_2;
        const auto V_MB   = (V_1 + V_2 -
            std::sqrt(V_diff * V_diff + 4 * this->delta_sq_)) / 2;
        const auto chi = std::log((V_MB - V_1) / this->delta_);

        std::ostringstream oss;
        oss << std::setw(18) << std::fixed << std::right << V_MB << ' ';
        oss << std::setw(16) << std::fixed << std::right << chi  << ' ';
        fmt += oss.str();

        return V_MB;
    }

    // for tests
    real_type delta() const noexcept {return delta_;}
    real_type dV1()   const noexcept {return dV1_;}
    real_type dV2()   const noexcept {return dV2_;}

    std::string const& name1() const noexcept {return name1_;}
    std::string const& name2() const noexcept {return name2_;}

    local_forcefield_type    const& local1()    const noexcept {return loc1_;}
    local_forcefield_type    const& local2()    const noexcept {return loc2_;}
    global_forcefield_type   const& global1()   const noexcept {return glo1_;}
    global_forcefield_type   const& global2()   const noexcept {return glo2_;}
    external_forcefield_type const& external1() const noexcept {return ext1_;}
    external_forcefield_type const& external2() const noexcept {return ext2_;}

    // -----------------------------------------------------------------------
    // calc_force/energy for each basin

    real_type calc_force_and_energy_basin1(system_type& sys) const
    {
        real_type energy = 0;
        energy += loc1_.calc_force_and_energy(sys);
        energy += glo1_.calc_force_and_energy(sys);
        energy += ext1_.calc_force_and_energy(sys);
        return energy;
    }
    real_type calc_force_and_energy_basin2(system_type& sys) const
    {
        real_type energy = 0;
        energy += loc2_.calc_force_and_energy(sys);
        energy += glo2_.calc_force_and_energy(sys);
        energy += ext2_.calc_force_and_energy(sys);
        return energy;
    }

    real_type calc_energy_basin1(const system_type& sys) const
    {
        return loc1_.calc_energy(sys) + glo1_.calc_energy(sys) +
               ext1_.calc_energy(sys);
    }
    real_type calc_energy_basin2(const system_type& sys) const
    {
        return loc2_.calc_energy(sys) + glo2_.calc_energy(sys) +
               ext2_.calc_energy(sys);
    }

  private:

    real_type dV1_;
    real_type dV2_;
    real_type delta_;
    real_type rdelta_;
    real_type delta_sq_;

    // this will be calculated in calc_energy(), that is marked as const.
    mutable real_type c2_over_c1_;

    std::string              name1_, name2_;
    local_forcefield_type    loc1_, loc2_;
    global_forcefield_type   glo1_, glo2_;
    external_forcefield_type ext1_, ext2_;

    // since it is used to contain temporary force from basin1,
    // it will be modified in calc_force() that is marked as const.
    // take care.
    mutable coordinate_container_type force_buffer_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class MultipleBasin2BasinUnit<SimulatorTraits<double, UnlimitedBoundary       >>;
extern template class MultipleBasin2BasinUnit<SimulatorTraits<float,  UnlimitedBoundary       >>;
extern template class MultipleBasin2BasinUnit<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class MultipleBasin2BasinUnit<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
