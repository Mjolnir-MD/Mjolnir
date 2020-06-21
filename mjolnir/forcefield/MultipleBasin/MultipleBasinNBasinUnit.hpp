#ifndef MJOLNIR_CORE_MULTIPLE_BASIN_N_BASIN_FORCE_FIELD_HPP
#define MJOLNIR_CORE_MULTIPLE_BASIN_N_BASIN_FORCE_FIELD_HPP
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/math/vector_util.hpp>
#include <mjolnir/util/string.hpp>
#include <algorithm>
#include <numeric>
#include <memory>

namespace mjolnir
{

// MultipleBasin forcefield with 4 or more basin.
// For 2 basin and 3 basin cases, see MultipleBasin2BasinForcefield.hpp and
// 3Basin one.
// If the number of basins are larger than 3, only 2 minimum-energy
// potentials are considered in the same way as the 2-basin case.
// It is the same strategy as the reference implementation.
//
// TODO: to speedup this forcefield...
// - add `calc_force_and_energy` member function to interactions
//   - it is technically easy. but it requires a huge effort.
//
template<typename traitsT>
class MultipleBasinNBasinForceField : public ForceFieldBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = ForceFieldBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using local_forcefield_type    = LocalForceField<traits_type>;
    using global_forcefield_type   = GlobalForceField<traits_type>;
    using external_forcefield_type = ExternalForceField<traits_type>;
    using forcefield_type          = std::tuple<
        local_forcefield_type, global_forcefield_type, external_forcefield_type>;
    using coordinate_container_type =
        typename system_type::coordinate_container_type;

  public:

    MultipleBasinNBasinForceField(
        std::vector<real_type>&& deltas, std::vector<real_type>&& dVs,
        std::vector<forcefield_type>&& basins, forcefield_type&& common)
        : deltas_(std::move(deltas)), dVs_(std::move(dVs)),
          basins_(std::move(basins)), common_(std::move(common))
    {}

    ~MultipleBasinNBasinForceField() override = default;
    MultipleBasinNBasinForceField(const MultipleBasinNBasinForceField&) = delete;
    MultipleBasinNBasinForceField(MultipleBasinNBasinForceField&&)      = default;
    MultipleBasinNBasinForceField& operator=(const MultipleBasinNBasinForceField&) = delete;
    MultipleBasinNBasinForceField& operator=(MultipleBasinNBasinForceField&&)      = default;

//     void initialize(system_type& sys) override
//     {
//         MJOLNIR_GET_DEFAULT_LOGGER();
//         MJOLNIR_LOG_FUNCTION();
//
//         this->force_buffer_.resize(sys.size(),
//                 math::make_coordinate<coordinate_type>(0, 0, 0));
//
//         // -------------------------------------------------------------------
//         MJOLNIR_LOG_INFO("initializing topology of the first basin");
//
//         loc1_         .write_topology(this->topol1_);
//         loc_common_   .write_topology(this->topol1_); // don't forget this
//         this->topol1_.construct_molecules();
//         loc1_.initialize(sys);
//         glo1_.initialize(sys, topol1_);
//         ext1_.initialize(sys);
//
//         // -------------------------------------------------------------------
//         MJOLNIR_LOG_INFO("initializing topology of the first basin");
//
//         loc2_         .write_topology(this->topol2_);
//         loc_common_   .write_topology(this->topol2_); // don't forget this
//         this->topol2_.construct_molecules();
//         loc2_.initialize(sys);
//         glo2_.initialize(sys, topol2_);
//         ext2_.initialize(sys);
//
//         // -------------------------------------------------------------------
//         // initialize the common part.
//
//         if(topol1_ == topol2_)
//         {
//             MJOLNIR_LOG_INFO("The topologies of two basins are the same.");
//             this->topologies_are_the_same_ = true;
//             loc_common_.initialize(sys);
//             glo_common_.initialize(sys, this->topol1_); // can use any one of them.
//             ext_common_.initialize(sys);
//         }
//         else
//         {
//             MJOLNIR_LOG_WARN("Two basins in MultipleBasin have different "
//                 "topologies. The common parts of global forcefield will only "
//                 "considers the topology defined in the common part of the "
//                 "local forcefield.");
//
//             topology_type topol_common_part_only;
//             topol_common_part_only.resize(sys.size());
//             loc_common_.write_topology(topol_common_part_only);
//
//             loc_common_.initialize(sys);
//             glo_common_.initialize(sys, this->topol_common_only_);
//             ext_common_.initialize(sys);
//         }
//         return;
//     }
//
//     void calc_force(system_type& sys) const noexcept override
//     {
//         // -------------------------------------------------------------------
//         // calc force of V_MB first.
//
//         this->calc_force_basin1(sys);
//         {
//             // save the current forces to force_buffer_.
//             // force_buffer is zero-cleared (at the end of this function),
//             // so the forces in the system will be zero-cleared after this.
//             using std::swap;
//             swap(this->force_buffer_, sys.forces());
//         }
//         this->calc_force_basin2(sys);
//
//         // TODO add calc_force_and_energy() to `Interaction`s.
//         const auto V_1    = this->calc_energy_basin1(sys) + this->dV1_;
//         const auto V_2    = this->calc_energy_basin2(sys) + this->dV2_;
//
//         // V_c does not affect to the diff because those are common.
//         // (V_1 + V_c) - (V_2 + V_c) = (V_1 - V_2)
//         const auto V_diff = V_1 - V_2;
//         const auto coef   = V_diff / std::sqrt(V_diff * V_diff + 4 * delta_sq_);
//
//         // F_V_MB = 1/2 [1 - Vdiff / sqrt(Vdiff^2 + 4 delta^2)] * F_V_1 +
//         //          1/2 [1 + Vdiff / sqrt(Vdiff^2 + 4 delta^2)] * F_V_2
//
//         const auto coef1 = real_type(0.5) * (real_type(1) - coef);
//         const auto coef2 = real_type(0.5) * (real_type(1) + coef);
//
//         // here, sys.forces has forces of basin2. force_buffer has forces of 1.
//         for(std::size_t i=0; i<sys.size(); ++i)
//         {
//             sys.force(i) *= coef2;
//             sys.force(i) += coef1 * force_buffer_[i];
//             force_buffer_[i] = math::make_coordinate<coordinate_type>(0, 0, 0);
//         }
//
//         // -------------------------------------------------------------------
//         // add forces from common part.
//         this->calc_force_common(sys);
//
//         return ;
//     }
//     real_type calc_energy(const system_type& sys) const noexcept override
//     {
//         // -------------------------------------------------------------------
//         // calc energy of V_MB
//         const auto V_1    = this->calc_energy_basin1(sys) + this->dV1_;
//         const auto V_2    = this->calc_energy_basin2(sys) + this->dV2_;
//         const auto V_diff = V_1 - V_2;
//
//         const auto V_MB =
//             (V_1 + V_2 - std::sqrt(V_diff * V_diff + 4 * delta_sq_)) / 2;
//
//         this->c2_over_c1_ = (V_MB - V_1) * this->rdelta_;
//
//         return V_MB + this->calc_energy_common(sys);
//     }
//
//     void update(const system_type& sys) override
//     {
//         // update parameters (e.g. temperature). TODO: topologies?
//
//         loc1_.update(sys);
//         glo1_.update(sys);
//         ext1_.update(sys);
//
//         loc2_.update(sys);
//         glo2_.update(sys);
//         ext2_.update(sys);
//
//         loc_common_.update(sys);
//         glo_common_.update(sys);
//         ext_common_.update(sys);
//         return;
//     }
//
//     // update margin of neighbor list
//     void update_margin(const real_type dmargin, const system_type& sys) override
//     {
//         loc1_.update_margin(dmargin, sys);
//         glo1_.update_margin(dmargin, sys);
//         ext1_.update_margin(dmargin, sys);
//
//         loc2_.update_margin(dmargin, sys);
//         glo2_.update_margin(dmargin, sys);
//         ext2_.update_margin(dmargin, sys);
//
//         loc_common_.update_margin(dmargin, sys);
//         glo_common_.update_margin(dmargin, sys);
//         ext_common_.update_margin(dmargin, sys);
//         return;
//     }
//
//     std::vector<std::string> list_energy_name() const override
//     {
//         using namespace mjolnir::literals::string_literals;
//         // Basin1{BondLength:Harmonic ...} Basin2{...} Common{...} Total chi
//
//         std::vector<std::string> basin1 = this->list_energy_name_basin1();
//         if(!basin1.empty())
//         {
//             basin1.front() = "Basin1{"_s + basin1.front();
//             basin1.back() += "}"_s;
//         }
//         std::vector<std::string> basin2 = this->list_energy_name_basin2();
//         if(!basin2.empty())
//         {
//             basin2.front() = "Basin2{"_s + basin2.front();
//             basin2.back() += "}"_s;
//         }
//         std::vector<std::string> common = this->list_energy_name_common();
//         if(!common.empty())
//         {
//             common.front() = "Common{"_s + common.front();
//             common.back() += "}"_s;
//         }
//
//         std::vector<std::string> retval;
//         retval.reserve(basin1.size() + basin2.size() + common.size() + 2);
//
//         std::copy(std::make_move_iterator(basin1.begin()),
//             std::make_move_iterator(basin1.end()), std::back_inserter(retval));
//         std::copy(std::make_move_iterator(basin2.begin()),
//             std::make_move_iterator(basin2.end()), std::back_inserter(retval));
//         std::copy(std::make_move_iterator(common.begin()),
//             std::make_move_iterator(common.end()), std::back_inserter(retval));
//
//         retval.push_back("MultipleBasinTotal");
//         retval.push_back("MultipleBasinChi");
//
//         return retval;
//     }
//     std::vector<real_type> dump_energy(const system_type& sys) const override
//     {
//         const auto es_1 = this->dump_energy_basin1(sys);
//         const auto es_2 = this->dump_energy_basin2(sys);
//         const auto es_c = this->dump_energy_common(sys);
//
//         const auto V_1 = std::accumulate(es_1.begin(), es_1.end(), real_type(0)) + this->dV1_;
//         const auto V_2 = std::accumulate(es_2.begin(), es_2.end(), real_type(0)) + this->dV2_;
//
//         const auto V_diff = V_1 - V_2;
//
//         const auto V_MB = (V_1 + V_2 -
//             std::sqrt(V_diff * V_diff + 4 * this->delta_sq_)) / 2;
//
//         // ( V1-V_MB   delta    ) (c1) = (0)
//         // (  delta  V2+dV-V_MB ) (c2)   (0)
//         // <=> (V1-V_MB) c1 + delta c2 = 0
//         // <=> (V_MB-V1) c1 = delta c2
//
//         const auto chi = std::log((V_MB - V_1) / this->delta_);
//
//         std::vector<real_type> retval;
//         retval.reserve(es_1.size() + es_2.size() + es_c.size() + 2);
//
//         std::copy(es_1.begin(), es_1.end(), std::back_inserter(retval));
//         std::copy(es_2.begin(), es_2.end(), std::back_inserter(retval));
//         std::copy(es_c.begin(), es_c.end(), std::back_inserter(retval));
//
//         retval.push_back(V_MB);
//         retval.push_back(chi);
//
//         return retval;
//     }
//
//     topology_type const& topology() const noexcept override
//     {
//         // currently, only DCDObserver uses this to write the number of chains.
//         if(this->topologies_are_the_same_)
//         {
//             return topol1_;
//         }
//         else if(this->c2_over_c1_ < real_type(1))
//         {
//             return topol1_;
//         }
//         else
//         {
//             return topol2_;
//         }
//     }
//
//     real_type delta() const noexcept {return delta_;}
//     real_type dV1()   const noexcept {return dV1_;}
//     real_type dV2()   const noexcept {return dV2_;}
//
//     // -----------------------------------------------------------------------
//     // calc_force/energy, dump/list_energy for each basin
//
//     void calc_force_basin1(system_type& sys) const
//     {
//         loc1_.calc_force(sys);
//         glo1_.calc_force(sys);
//         ext1_.calc_force(sys);
//         return;
//     }
//     void calc_force_basin2(system_type& sys) const
//     {
//         loc2_.calc_force(sys);
//         glo2_.calc_force(sys);
//         ext2_.calc_force(sys);
//         return;
//     }
//     void calc_force_common(system_type& sys) const
//     {
//         loc_common_.calc_force(sys);
//         glo_common_.calc_force(sys);
//         ext_common_.calc_force(sys);
//         return;
//     }
//
//     real_type calc_energy_basin1(const system_type& sys) const
//     {
//         return loc1_.calc_energy(sys) + glo1_.calc_energy(sys) +
//                ext1_.calc_energy(sys);
//     }
//     real_type calc_energy_basin2(const system_type& sys) const
//     {
//         return loc2_.calc_energy(sys) + glo2_.calc_energy(sys) +
//                ext2_.calc_energy(sys);
//     }
//     real_type calc_energy_common(const system_type& sys) const
//     {
//         return loc_common_.calc_energy(sys) + glo_common_.calc_energy(sys) +
//                ext_common_.calc_energy(sys);
//     }
//
//     std::vector<std::string> list_energy_name_basin1() const
//     {
//         auto retval = loc1_.list_energy();
//         auto glo    = glo1_.list_energy();
//         auto ext    = ext1_.list_energy();
//
//         retval.reserve(retval.size() + glo.size() + ext.size());
//
//         std::copy(std::make_move_iterator(glo.begin()),
//                   std::make_move_iterator(glo.end()),
//                   std::back_inserter(retval));
//
//         std::copy(std::make_move_iterator(ext.begin()),
//                   std::make_move_iterator(ext.end()),
//                   std::back_inserter(retval));
//         return retval;
//     }
//     std::vector<std::string> list_energy_name_basin2() const
//     {
//         auto retval = loc2_.list_energy();
//         auto glo    = glo2_.list_energy();
//         auto ext    = ext2_.list_energy();
//
//         retval.reserve(retval.size() + glo.size() + ext.size());
//
//         std::copy(std::make_move_iterator(glo.begin()),
//                   std::make_move_iterator(glo.end()),
//                   std::back_inserter(retval));
//
//         std::copy(std::make_move_iterator(ext.begin()),
//                   std::make_move_iterator(ext.end()),
//                   std::back_inserter(retval));
//         return retval;
//     }
//     std::vector<std::string> list_energy_name_common() const
//     {
//         auto retval = loc_common_.list_energy();
//         auto glo    = glo_common_.list_energy();
//         auto ext    = ext_common_.list_energy();
//
//         retval.reserve(retval.size() + glo.size() + ext.size());
//
//         std::copy(std::make_move_iterator(glo.begin()),
//                   std::make_move_iterator(glo.end()),
//                   std::back_inserter(retval));
//
//         std::copy(std::make_move_iterator(ext.begin()),
//                   std::make_move_iterator(ext.end()),
//                   std::back_inserter(retval));
//         return retval;
//     }
//
//     std::vector<real_type> dump_energy_basin1(const system_type& sys) const
//     {
//         auto retval = loc1_.dump_energy(sys);
//         auto glo    = glo1_.dump_energy(sys);
//         auto ext    = ext1_.dump_energy(sys);
//
//         retval.reserve(retval.size() + glo.size() + ext.size());
//
//         std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
//         std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
//         return retval;
//     }
//     std::vector<real_type> dump_energy_basin2(const system_type& sys) const
//     {
//         auto retval = loc2_.dump_energy(sys);
//         auto glo    = glo2_.dump_energy(sys);
//         auto ext    = ext2_.dump_energy(sys);
//
//         retval.reserve(retval.size() + glo.size() + ext.size());
//
//         std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
//         std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
//         return retval;
//     }
//     std::vector<real_type> dump_energy_common(const system_type& sys) const
//     {
//         auto retval = loc_common_.dump_energy(sys);
//         auto glo    = glo_common_.dump_energy(sys);
//         auto ext    = ext_common_.dump_energy(sys);
//
//         retval.reserve(retval.size() + glo.size() + ext.size());
//
//         std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
//         std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
//         return retval;
//     }

  private:

    std::vector<real_type> dVs_;
    std::vector<real_type> deltas_;
    std::vector<real_type> delta_sqs_;

    bool topologies_are_the_same_;

    std::vector<topology_type>   topols_;
    std::vector<forcefield_type> basins_;

    // since it is used to contain temporary force from basin1,
    // it will be modified in calc_force() that is marked as const.
    // take care.
    mutable std::vector<coordinate_container_type> force_buffers_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class MultipleBasinNBasinForceField<SimulatorTraits<double, UnlimitedBoundary       >>;
extern template class MultipleBasinNBasinForceField<SimulatorTraits<float,  UnlimitedBoundary       >>;
extern template class MultipleBasinNBasinForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class MultipleBasinNBasinForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
