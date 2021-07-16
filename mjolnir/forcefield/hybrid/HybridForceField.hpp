#ifndef MJOLNIR_FORCEFIELD_HYBRID_FORCE_FIELD_HPP
#define MJOLNIR_FORCEFIELD_HYBRID_FORCE_FIELD_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/Topology.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/string.hpp>
#include <algorithm>
#include <numeric>
#include <memory>

namespace mjolnir
{

// HybridForceField
//
// U = lambda U_1 + (1 - lambra) U_2
//
template<typename traitsT>
class HybridForceField : public ForceFieldBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = ForceFieldBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using coordinate_container_type = typename system_type::coordinate_container_type;
    using topology_type   = Topology;
    using forcefield_type = std::unique_ptr<base_type>;
    using constraint_forcefield_type  = ConstraintForceField<traits_type>;

  public:

    HybridForceField(const real_type lambda,
                     forcefield_type ff_1, forcefield_type ff_2)
        : lambda_(lambda), ff_1_(std::move(ff_1)), ff_2_(std::move(ff_2))
    {}

    ~HybridForceField() override = default;
    HybridForceField(const HybridForceField&) = delete;
    HybridForceField(HybridForceField&&)      = default;
    HybridForceField& operator=(const HybridForceField&) = delete;
    HybridForceField& operator=(HybridForceField&&)      = default;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();

        this->ff_1_->initialize(sys);
        this->ff_2_->initialize(sys);

        this->force_buffer_.resize(sys.size(),
                math::make_coordinate<coordinate_type>(0, 0, 0));

        if( ! ff_1_->constraint().empty() || ! ff_2_->constraint().empty())
        {
            MJOLNIR_LOG_ERROR("hybrid constraint cannot mix 2 different constraints");
        }
        return;
    }

    void calc_force(system_type& sys) const noexcept override
    {
        using std::swap;

        sys.preprocess_forces();
        this->ff_1_->calc_force(sys);
        sys.postprocess_forces();

        swap(this->force_buffer_, sys.forces());

        sys.preprocess_forces();
        this->ff_2_->calc_force(sys);
        sys.postprocess_forces();

        const real_type one_minus_lambda = real_type(1) - this->lambda_;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) *= one_minus_lambda;
            sys.force(i) += lambda_ * force_buffer_[i];
            force_buffer_[i] = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        return ;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        return lambda_ * ff_1_->calc_energy(sys) +
               (real_type(1) - lambda_) * ff_2_->calc_energy(sys);
    }

    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        using std::swap;

        sys.preprocess_forces();
        const auto V1 = this->ff_1_->calc_force_and_energy(sys);
        sys.postprocess_forces();

        swap(this->force_buffer_, sys.forces());

        sys.preprocess_forces();
        const auto V2 = this->ff_2_->calc_force_and_energy(sys);
        sys.postprocess_forces();

        const real_type one_minus_lambda = real_type(1) - this->lambda_;
        for(std::size_t i=0; i<sys.size(); ++i)
        {
            sys.force(i) *= one_minus_lambda;
            sys.force(i) += lambda_ * force_buffer_[i];
            force_buffer_[i] = math::make_coordinate<coordinate_type>(0, 0, 0);
        }
        return lambda_ * V1 + one_minus_lambda * V2;
    }

    void update(const system_type& sys) override
    {
        this->ff_1_->update(sys);
        this->ff_2_->update(sys);
        return;
    }

    // update margin of neighbor list
    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->ff_1_->reduce_margin(dmargin, sys);
        this->ff_2_->reduce_margin(dmargin, sys);
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->ff_1_->scale_margin(scale, sys);
        this->ff_2_->scale_margin(scale, sys);
        return;
    }

    constraint_forcefield_type const& constraint() const noexcept override
    {
        return ff_1_->constraint();
    }

    // This is actually not correct, but it is used only by DCDObserver to write
    // number of chains, so it should not be problematic (but need to be refactored)
    topology_type const& topology() const noexcept override {return ff_1_->topology();}

    real_type lambda() const noexcept {return lambda_;}

    forcefield_type const& ff1() const noexcept {return ff_1_;}
    forcefield_type const& ff2() const noexcept {return ff_2_;}

    // -----------------------------------------------------------------------
    // energy output format

    void format_energy_name(std::string& fmt) const override
    {
        using namespace mjolnir::literals::string_literals;

        fmt += "FF 1{"_s;
        ff_1_->format_energy_name(fmt);
        if(!fmt.empty() && fmt.back() == ' ') {fmt.pop_back();}
        fmt += "} "_s;

        fmt += "FF 2{"_s;
        ff_2_->format_energy_name(fmt);
        if(!fmt.empty() && fmt.back() == ' ') {fmt.pop_back();}
        fmt += "} "_s;
        return;
    }

    real_type format_energy(const system_type& sys, std::string& fmt) const override
    {
        real_type total = 0.0;

        fmt += "     "_s; // FF 1{
        total += ff_1_->format_energy(sys, fmt);
        if(!fmt.empty() && fmt.back() == ' ') {fmt.pop_back();}
        fmt += "  "_s; // }

        fmt += "     "_s; // FF 2{
        total += ff_2_->format_energy(sys, fmt);
        if(!fmt.empty() && fmt.back() == ' ') {fmt.pop_back();}
        fmt += "  "_s; // }
        return total;
    }

  private:

    real_type lambda_;
    mutable coordinate_container_type force_buffer_; // this will be used in calc_force
    forcefield_type ff_1_;
    forcefield_type ff_2_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class HybridForceField<SimulatorTraits<double, UnlimitedBoundary       >>;
extern template class HybridForceField<SimulatorTraits<float,  UnlimitedBoundary       >>;
extern template class HybridForceField<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class HybridForceField<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
#endif

} // mjolnir
#endif// MJOLNIR_CORE_MULTIPLE_BASIN_FORCE_FIELD_HPP
