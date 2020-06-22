#ifndef MJOLNIR_CORE_MULTIPLE_BASIN_3_BASIN_FORCE_FIELD_HPP
#define MJOLNIR_CORE_MULTIPLE_BASIN_3_BASIN_FORCE_FIELD_HPP
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

// 3-basin MultipleBasin forcefield.
//
// TODO: to speedup this forcefield...
// - add `calc_force_and_energy` member function to interactions
//   - it is technically easy. but it requires a huge effort.
//
// V_MB is a solution of the following equation.
//
// (V1 + dV1  delta12  delta13 ) (c1)        (c1)
// (delta21   V2 + dV2 delta23 ) (c2) = V_MB (c2)
// (delta31   delta32  V3 + dV3) (c3)        (c3)
//
// Here, after replacing Vi + dVi = Ui, the value of V_MB can be written as:
//
// 0 = -V_MB^3 + (U1+U2+U3) V_MB^2 -
//     (U1U2 + U2U3 + U3U1 + delta12^2 + delta23^2 + delta31^2) V_MB +
//     U1U2U3 + 2*delta12*delta23*delta31 +
//     U1*delta23^2 + U2*delta31^2 + U3*delta12^2
//
template<typename traitsT>
class MultipleBasin3BasinUnit final: public MultipleBasinUnitBase<traitsT>
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

    MultipleBasin3BasinUnit(
        const std::string& name1, const std::string& name2, const std::string& name3,
        const real_type  delta12, const real_type  delta23, const real_type  delta31,
        const real_type  dV1,     const real_type  dV2,     const real_type  dV3,
        forcefield_type  basin1,  forcefield_type  basin2,  forcefield_type  basin3)
        : dV1_(dV1),  dV2_(dV2),  dV3_(dV3),
          delta12_(delta12), delta23_(delta23), delta31_(delta31),
          delta12_sq_(delta12 * delta12), delta23_sq_(delta23 * delta23),
          delta31_sq_(delta31 * delta31), c1_(1), c2_(1), c3_(1),
          name1_(name1), name2_(name2), name3_(name3),
          basin1_(std::move(basin1)), basin2_(std::move(basin2)),
          basin3_(std::move(basin3))
    {}

    ~MultipleBasin3BasinUnit() override = default;
    MultipleBasin3BasinUnit(const MultipleBasin3BasinUnit&) = delete;
    MultipleBasin3BasinUnit(MultipleBasin3BasinUnit&&)      = default;
    MultipleBasin3BasinUnit& operator=(const MultipleBasin3BasinUnit&) = delete;
    MultipleBasin3BasinUnit& operator=(MultipleBasin3BasinUnit&&)      = default;

    void write_topology(const system_type& sys, topology_type& topol) const override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        MJOLNIR_LOG_INFO("checking topologies are the same");

        Topology topol1(sys.size()), topol2(sys.size()), topol3(sys.size());

        std::get<0>(basin1_).write_topology(topol1);
        std::get<0>(basin2_).write_topology(topol2);
        std::get<0>(basin3_).write_topology(topol3);
        topol1.construct_molecules();
        topol2.construct_molecules();
        topol3.construct_molecules();
        if(topol1 != topol2)
        {
            MJOLNIR_LOG_ERROR("topologies of 2 basins (", name1_, " and ",
                              name2_, ") are different from each other.");
            MJOLNIR_LOG_ERROR("MultipleBasin does not support such a case.");
            throw std::runtime_error("mjolnir::MultipleBasin2BasinUnit: "
                    "Topologies of 2 basins shouold be the same");
        }
        if(topol2 != topol3)
        {
            MJOLNIR_LOG_ERROR("topologies of 2 basins (", name2_, " and ",
                              name3_, ") are different from each other.");
            MJOLNIR_LOG_ERROR("MultipleBasin does not support such a case.");
            throw std::runtime_error("mjolnir::MultipleBasin2BasinUnit: "
                    "Topologies of 2 basins shouold be the same");
        }

        MJOLNIR_LOG_INFO("writing topology");
        std::get<0>(basin1_).write_topology(topol);
        return;
    }

    void initialize(const system_type& sys, const topology_type& topol) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->force_buffer1_.resize(sys.size(),
                math::make_coordinate<coordinate_type>(0, 0, 0));
        this->force_buffer2_.resize(sys.size(),
                math::make_coordinate<coordinate_type>(0, 0, 0));

        std::get<0>(basin1_).initialize(sys);
        std::get<1>(basin1_).initialize(sys, topol);
        std::get<2>(basin1_).initialize(sys);

        std::get<0>(basin2_).initialize(sys);
        std::get<1>(basin2_).initialize(sys, topol);
        std::get<2>(basin2_).initialize(sys);

        std::get<0>(basin3_).initialize(sys);
        std::get<1>(basin3_).initialize(sys, topol);
        std::get<2>(basin3_).initialize(sys);

        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        using std::swap;
        // -------------------------------------------------------------------
        // calc force of V_MB first.

        // save the current forces to force_buffer_.
        // force_buffer is zero-cleared (at the end of this function),
        // so the forces in the system will be zero-cleared after this.
        this->calc_force_basin1(sys);
        swap(this->force_buffer1_, sys.forces());

        this->calc_force_basin2(sys);
        swap(this->force_buffer2_, sys.forces());

        this->calc_force_basin3(sys);

        // TODO add calc_force_and_energy() to `Interaction`s.
        const auto V_1 = this->calc_energy_basin1(sys) + this->dV1_;
        const auto V_2 = this->calc_energy_basin2(sys) + this->dV2_;
        const auto V_3 = this->calc_energy_basin3(sys) + this->dV3_;

        const auto V_MB = this->calc_V_MB(V_1, V_2, V_3);

        const auto V_diff_1 = V_1 - V_MB;
        const auto V_diff_2 = V_2 - V_MB;
        const auto V_diff_3 = V_3 - V_MB;
        const auto denom = real_type(1) /
            (V_diff_1 * V_diff_2 + V_diff_2 * V_diff_3 + V_diff_3 * V_diff_1 -
             delta12_sq_ - delta23_sq_ - delta31_sq_);

        const auto coef1 = (V_diff_2 * V_diff_3 - delta23_sq_) * denom;
        const auto coef2 = (V_diff_3 * V_diff_1 - delta31_sq_) * denom;
        const auto coef3 = (V_diff_1 * V_diff_2 - delta12_sq_) * denom;

        // here, sys.forces has forces of basin2. force_buffer has forces of 1.
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) *= coef3;
            sys.force(i) += coef1 * force_buffer1_[i] +
                            coef2 * force_buffer2_[i];
            force_buffer1_[i] = math::make_coordinate<coordinate_type>(0, 0, 0);
            force_buffer2_[i] = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        return ;
    }
    real_type calc_energy(const system_type& sys) const noexcept override
    {
        // -------------------------------------------------------------------
        // calc energy of V_MB
        const auto V_1  = this->calc_energy_basin1(sys) + this->dV1_;
        const auto V_2  = this->calc_energy_basin2(sys) + this->dV2_;
        const auto V_3  = this->calc_energy_basin3(sys) + this->dV3_;
        const auto V_MB = this->calc_V_MB(V_1, V_2, V_3);

        this->c1_ = 1; // always 1, baseline.
        this->c2_ = (delta12_ * (V_MB - V_3) + delta23_ * delta31_) /
                    ((V_MB - V_2) * (V_MB - V_3) - delta23_sq_);
        this->c3_ = (delta31_ * (V_MB - V_2) + delta12_ * delta23_) /
                    ((V_MB - V_2) * (V_MB - V_3) - delta23_sq_);

        return V_MB;
    }

    void update(const system_type& sys, const topology_type& topol) override
    {
        // update parameters (e.g. temperature). TODO: topologies?

        std::get<0>(basin1_).update(sys);
        std::get<1>(basin1_).update(sys, topol);
        std::get<2>(basin1_).update(sys);

        std::get<0>(basin2_).update(sys);
        std::get<1>(basin2_).update(sys, topol);
        std::get<2>(basin2_).update(sys);

        std::get<0>(basin3_).update(sys);
        std::get<1>(basin3_).update(sys, topol);
        std::get<2>(basin3_).update(sys);
        return;
    }

    // update margin of neighbor list
    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        std::get<0>(basin1_).reduce_margin(dmargin, sys);
        std::get<1>(basin1_).reduce_margin(dmargin, sys);
        std::get<2>(basin1_).reduce_margin(dmargin, sys);

        std::get<0>(basin2_).reduce_margin(dmargin, sys);
        std::get<1>(basin2_).reduce_margin(dmargin, sys);
        std::get<2>(basin2_).reduce_margin(dmargin, sys);

        std::get<0>(basin3_).reduce_margin(dmargin, sys);
        std::get<1>(basin3_).reduce_margin(dmargin, sys);
        std::get<2>(basin3_).reduce_margin(dmargin, sys);
        return;
    }
    void scale_margin(const real_type dmargin, const system_type& sys) override
    {
        std::get<0>(basin1_).scale_margin(dmargin, sys);
        std::get<1>(basin1_).scale_margin(dmargin, sys);
        std::get<2>(basin1_).scale_margin(dmargin, sys);

        std::get<0>(basin2_).scale_margin(dmargin, sys);
        std::get<1>(basin2_).scale_margin(dmargin, sys);
        std::get<2>(basin2_).scale_margin(dmargin, sys);

        std::get<0>(basin3_).scale_margin(dmargin, sys);
        std::get<1>(basin3_).scale_margin(dmargin, sys);
        std::get<2>(basin3_).scale_margin(dmargin, sys);
        return;
    }

    std::vector<std::string> list_energy_name() const override
    {
        using namespace mjolnir::literals::string_literals;
        // Basin1{BondLength:Harmonic ...} Basin2{...} Basin3{...} Total c1 c2 c3

        std::vector<std::string> basin1 = this->list_energy_name_basin1();
        if(!basin1.empty())
        {
            basin1.front() = name1_ + "{"_s + basin1.front();
            basin1.back() += "}"_s;
        }
        std::vector<std::string> basin2 = this->list_energy_name_basin2();
        if(!basin2.empty())
        {
            basin2.front() = name2_ + "{"_s + basin2.front();
            basin2.back() += "}"_s;
        }
        std::vector<std::string> basin3 = this->list_energy_name_basin3();
        if(!basin3.empty())
        {
            basin3.front() = name3_ + "{"_s + basin3.front();
            basin3.back() += "}"_s;
        }

        std::vector<std::string> retval;
        retval.reserve(basin1.size() + basin2.size() + basin3.size() + 4);

        std::copy(std::make_move_iterator(basin1.begin()),
            std::make_move_iterator(basin1.end()), std::back_inserter(retval));
        std::copy(std::make_move_iterator(basin2.begin()),
            std::make_move_iterator(basin2.end()), std::back_inserter(retval));
        std::copy(std::make_move_iterator(basin3.begin()),
            std::make_move_iterator(basin3.end()), std::back_inserter(retval));

        retval.push_back("MultipleBasinTotal");
        retval.push_back("MultipleBasin_c1");
        retval.push_back("MultipleBasin_c2");
        retval.push_back("MultipleBasin_c3");

        return retval;
    }
    std::vector<real_type> dump_energy(const system_type& sys) const override
    {
        const auto es_1 = this->dump_energy_basin1(sys);
        const auto es_2 = this->dump_energy_basin2(sys);
        const auto es_3 = this->dump_energy_basin3(sys);

        const auto V_1 = std::accumulate(es_1.begin(), es_1.end(), real_type(0)) + this->dV1_;
        const auto V_2 = std::accumulate(es_2.begin(), es_2.end(), real_type(0)) + this->dV2_;
        const auto V_3 = std::accumulate(es_3.begin(), es_3.end(), real_type(0)) + this->dV3_;

        const auto V_MB = this->calc_V_MB(V_1, V_2, V_3);

        this->c1_ = 1; // always 1, baseline.
        this->c2_ = (delta12_ * (V_MB - V_3) + delta23_ * delta31_) /
                    ((V_MB - V_2) * (V_MB - V_3) - delta23_sq_);
        this->c3_ = (delta31_ * (V_MB - V_2) + delta12_ * delta23_) /
                    ((V_MB - V_2) * (V_MB - V_3) - delta23_sq_);

        std::vector<real_type> retval;
        retval.reserve(es_1.size() + es_2.size() + es_3.size() + 4);

        std::copy(es_1.begin(), es_1.end(), std::back_inserter(retval));
        std::copy(es_2.begin(), es_2.end(), std::back_inserter(retval));
        std::copy(es_3.begin(), es_3.end(), std::back_inserter(retval));

        retval.push_back(V_MB);
        retval.push_back(this->c1_);
        retval.push_back(this->c2_);
        retval.push_back(this->c3_);
        return retval;
    }

    // for tests
    real_type delta12() const noexcept {return delta12_;}
    real_type delta23() const noexcept {return delta23_;}
    real_type delta31() const noexcept {return delta31_;}
    real_type dV1()     const noexcept {return dV1_;}
    real_type dV2()     const noexcept {return dV2_;}
    real_type dV3()     const noexcept {return dV3_;}

    std::string const& name1() const noexcept {return name1_;}
    std::string const& name2() const noexcept {return name2_;}
    std::string const& name3() const noexcept {return name3_;}

    local_forcefield_type    const& local1()    const noexcept {return std::get<0>(basin1_);}
    local_forcefield_type    const& local2()    const noexcept {return std::get<0>(basin2_);}
    local_forcefield_type    const& local3()    const noexcept {return std::get<0>(basin3_);}
    global_forcefield_type   const& global1()   const noexcept {return std::get<1>(basin1_);}
    global_forcefield_type   const& global2()   const noexcept {return std::get<1>(basin2_);}
    global_forcefield_type   const& global3()   const noexcept {return std::get<1>(basin3_);}
    external_forcefield_type const& external1() const noexcept {return std::get<2>(basin1_);}
    external_forcefield_type const& external2() const noexcept {return std::get<2>(basin2_);}
    external_forcefield_type const& external3() const noexcept {return std::get<2>(basin3_);}

    // -----------------------------------------------------------------------
    // calc_force/energy, dump/list_energy for each basin

    void calc_force_basin1(system_type& sys) const
    {
        std::get<0>(basin1_).calc_force(sys);
        std::get<1>(basin1_).calc_force(sys);
        std::get<2>(basin1_).calc_force(sys);
        return;
    }
    void calc_force_basin2(system_type& sys) const
    {
        std::get<0>(basin2_).calc_force(sys);
        std::get<1>(basin2_).calc_force(sys);
        std::get<2>(basin2_).calc_force(sys);
        return;
    }
    void calc_force_basin3(system_type& sys) const
    {
        std::get<0>(basin3_).calc_force(sys);
        std::get<1>(basin3_).calc_force(sys);
        std::get<2>(basin3_).calc_force(sys);
        return;
    }

    real_type calc_energy_basin1(const system_type& sys) const
    {
        return std::get<0>(basin1_).calc_energy(sys) +
               std::get<1>(basin1_).calc_energy(sys) +
               std::get<2>(basin1_).calc_energy(sys);
    }
    real_type calc_energy_basin2(const system_type& sys) const
    {
        return std::get<0>(basin2_).calc_energy(sys) +
               std::get<1>(basin2_).calc_energy(sys) +
               std::get<2>(basin2_).calc_energy(sys);
    }
    real_type calc_energy_basin3(const system_type& sys) const
    {
        return std::get<0>(basin3_).calc_energy(sys) +
               std::get<1>(basin3_).calc_energy(sys) +
               std::get<2>(basin3_).calc_energy(sys);
    }

    std::vector<std::string> list_energy_name_basin1() const
    {
        auto retval = std::get<0>(basin1_).list_energy();
        auto glo    = std::get<1>(basin1_).list_energy();
        auto ext    = std::get<2>(basin1_).list_energy();

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(std::make_move_iterator(glo.begin()),
                  std::make_move_iterator(glo.end()),
                  std::back_inserter(retval));
        std::copy(std::make_move_iterator(ext.begin()),
                  std::make_move_iterator(ext.end()),
                  std::back_inserter(retval));
        return retval;
    }
    std::vector<std::string> list_energy_name_basin2() const
    {
        auto retval = std::get<0>(basin2_).list_energy();
        auto glo    = std::get<1>(basin2_).list_energy();
        auto ext    = std::get<2>(basin2_).list_energy();

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(std::make_move_iterator(glo.begin()),
                  std::make_move_iterator(glo.end()),
                  std::back_inserter(retval));
        std::copy(std::make_move_iterator(ext.begin()),
                  std::make_move_iterator(ext.end()),
                  std::back_inserter(retval));
        return retval;
    }
    std::vector<std::string> list_energy_name_basin3() const
    {
        auto retval = std::get<0>(basin3_).list_energy();
        auto glo    = std::get<1>(basin3_).list_energy();
        auto ext    = std::get<2>(basin3_).list_energy();

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(std::make_move_iterator(glo.begin()),
                  std::make_move_iterator(glo.end()),
                  std::back_inserter(retval));
        std::copy(std::make_move_iterator(ext.begin()),
                  std::make_move_iterator(ext.end()),
                  std::back_inserter(retval));
        return retval;
    }

    std::vector<real_type> dump_energy_basin1(const system_type& sys) const
    {
        auto retval = std::get<0>(basin1_).dump_energy(sys);
        auto glo    = std::get<1>(basin1_).dump_energy(sys);
        auto ext    = std::get<2>(basin1_).dump_energy(sys);

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
        std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
        return retval;
    }
    std::vector<real_type> dump_energy_basin2(const system_type& sys) const
    {
        auto retval = std::get<0>(basin2_).dump_energy(sys);
        auto glo    = std::get<1>(basin2_).dump_energy(sys);
        auto ext    = std::get<2>(basin2_).dump_energy(sys);

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
        std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
        return retval;
    }
    std::vector<real_type> dump_energy_basin3(const system_type& sys) const
    {
        auto retval = std::get<0>(basin3_).dump_energy(sys);
        auto glo    = std::get<1>(basin3_).dump_energy(sys);
        auto ext    = std::get<2>(basin3_).dump_energy(sys);

        retval.reserve(retval.size() + glo.size() + ext.size());

        std::copy(glo.begin(), glo.end(), std::back_inserter(retval));
        std::copy(ext.begin(), ext.end(), std::back_inserter(retval));
        return retval;
    }

  private:

    real_type calc_V_MB(const real_type V_1, const real_type V_2,
                        const real_type V_3) const noexcept
    {
        const auto b = V_1 + V_2 + V_3;
        const auto c = 0.5 * (V_1*V_1 + V_2*V_2 + V_3*V_3 - b * b) +
                      delta12_sq_ + delta23_sq_ + delta31_sq_;
        const auto d = V_1 * V_2 * V_3 + 2 * delta12_ * delta23_ * delta31_ -
                       V_1 * delta23_sq_ - V_2 * delta31_sq_ - V_3 * delta12_sq_;

        return this->min_solution_of_cubic_equation(-b, -c, -d); // a == -1
    }

    // minimum solution of 0 = x^3 + bx^2 + cx + d.
    static real_type min_solution_of_cubic_equation(// assuming a == 1
            const real_type b, const real_type c, const real_type d) noexcept
    {
        constexpr real_type one_over_3  = real_type(1) / real_type(3);
        constexpr real_type one_over_27 = real_type(1) / real_type(27);

        constexpr real_type one_over_2  = real_type(1) / real_type(2);
        constexpr real_type one_over_4  = real_type(1) / real_type(4);

        constexpr real_type abs_tol = math::abs_tolerance<real_type>();
        constexpr real_type pi      = math::constants<real_type>::pi();

        const real_type bsq = b * b;

        const real_type p = c - bsq * one_over_3;
        const real_type q = d + b * (2 * bsq - 9 * c) * one_over_27;

        const real_type discriminant = one_over_4  * q * q +
                                       one_over_27 * p * p * p;

        real_type root_min = std::numeric_limits<real_type>::quiet_NaN();
        if(0 <= discriminant)
        {
            const auto sqrtD = std::sqrt(discriminant);
            const auto s11 = -one_over_2 * q + sqrtD;
            const auto s22 = -one_over_2 * q - sqrtD;

            const auto s1 = std::copysign(std::cbrt(std::abs(s11)), s11);
            const auto s2 = std::copysign(std::cbrt(std::abs(s22)), s22);

            root_min = s1 + s2;

            if(discriminant < abs_tol) // D == 0
            {
                root_min = std::min(root_min, -s1);
            }
        }
        else // 3 real roots, 2 are the same
        {
            const auto pabs = std::abs(p);
            const auto s1   = std::sqrt(pabs * pabs * pabs * one_over_27);
            const auto s2   = std::sqrt(pabs * one_over_3);
            const auto phi  = std::acos(-one_over_2 * q / s1);

            root_min = std::min(2 * s2 * std::cos(phi * one_over_3), std::min(
                        -2 * s2 * std::cos((phi + pi) * one_over_3),
                        -2 * s2 * std::cos((phi - pi) * one_over_3)));
        }
        return root_min - b * one_over_3;
    }

  private:

    real_type dV1_;
    real_type dV2_;
    real_type dV3_;
    real_type delta12_;
    real_type delta23_;
    real_type delta31_;
    real_type delta12_sq_;
    real_type delta23_sq_;
    real_type delta31_sq_;

    // those are calculated in calc_energy().
    mutable real_type c1_; // always 1.0, as a baseline (exists for consistency)
    mutable real_type c2_; // relative to c1
    mutable real_type c3_; // relative to c1

    std::string     name1_,  name2_,  name3_;
    forcefield_type basin1_, basin2_, basin3_;

    // XXX since it is used to contain temporary force from basin1/2,
    // XXX it will be modified in calc_force() that is marked as const.
    // XXX take care.
    mutable coordinate_container_type force_buffer1_;
    mutable coordinate_container_type force_buffer2_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class MultipleBasin3BasinUnit<SimulatorTraits<double, UnlimitedBoundary       >>;
extern template class MultipleBasin3BasinUnit<SimulatorTraits<float,  UnlimitedBoundary       >>;
extern template class MultipleBasin3BasinUnit<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class MultipleBasin3BasinUnit<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
