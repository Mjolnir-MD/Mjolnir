#ifndef MJOLNIR_INTERACTION_BOND_ANGLE_INTERACTION_HPP
#define MJOLNIR_INTERACTION_BOND_ANGLE_INTERACTION_HPP
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/math/constants.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <cmath>

namespace mjolnir
{

// calculate energy and force of bond angle interaction.
// The implementation is based on the book written by A. Ueda (2003) Shokabo.
template<typename traitsT, typename potentialT>
class BondAngleInteraction final : public LocalInteractionBase<traitsT>
{
  public:
    using traits_type          = traitsT;
    using potential_type       = potentialT;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using indices_type         = std::array<std::size_t, 3>;
    using potential_index_pair = std::pair<indices_type, potentialT>;
    using container_type       = std::vector<potential_index_pair>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

  public:

    BondAngleInteraction(const connection_kind_type kind,
                         const container_type& pot)
        : kind_(kind), potentials_(pot)
    {}
    BondAngleInteraction(const connection_kind_type kind,
                         container_type&& pot)
        : kind_(kind), potentials_(std::move(pot))
    {}
    BondAngleInteraction(const BondAngleInteraction&) = default;
    BondAngleInteraction(BondAngleInteraction&&)      = default;
    BondAngleInteraction& operator=(const BondAngleInteraction&) = default;
    BondAngleInteraction& operator=(BondAngleInteraction&&)      = default;
    ~BondAngleInteraction() override {}

    real_type calc_energy(const system_type&) const noexcept override;

    void calc_force(system_type& sys) const noexcept override
    {
        this->template calc_force_energy_virial_impl<false, false>(sys);
        return ;
    }
    void calc_force_and_virial(system_type& sys) const noexcept override
    {
        this->template calc_force_energy_virial_impl<false, true>(sys);
        return ;
    }
    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        return this->template calc_force_energy_virial_impl<true, true>(sys);
    }

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of angles = ", potentials_.size());
        for(auto& potential : this->potentials_)
        {
            potential.second.initialize(sys);
        }
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials_)
        {
            item.second.update(sys);
        }
    }

    // do nothing. this is used to reduce margin of neighbor list, and added
    // to this class for the consistency.
    void reduce_margin(const real_type, const system_type&) override {return;}
    void  scale_margin(const real_type, const system_type&) override {return;}

    std::string name() const override
    {return "BondAngle:"_s + potential_type::name();}

    void write_topology(topology_type&) const override;

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new BondAngleInteraction(kind_, container_type(potentials_));
    }

  private:

    template<bool NeedEnergy, bool NeedVirial>
    real_type calc_force_energy_virial_impl(system_type& sys) const noexcept;

  private:
    connection_kind_type kind_;
    container_type potentials_;

#ifdef MJOLNIR_WITH_OPENMP
    // OpenMP implementation uses its own implementation to run it in parallel.
    // So this implementation should not be instanciated with OpenMP Traits.
    static_assert(!is_openmp_simulator_traits<traits_type>::value,
                  "this is the default implementation, not for OpenMP");
#endif
};

template<typename traitsT, typename potentialT>
template<bool NeedEnergy, bool NeedVirial>
typename BondAngleInteraction<traitsT, potentialT>::real_type
BondAngleInteraction<traitsT, potentialT>::calc_force_energy_virial_impl(
        system_type& sys) const noexcept
{
    constexpr auto abs_tol = math::abs_tolerance<real_type>();

    real_type energy = 0;
    for(const auto& idxp : this->potentials_)
    {
        const auto& idxs = idxp.first;
        const auto& pot  = idxp.second;

        const auto& p0 = sys.position(idxs[0]);
        const auto& p1 = sys.position(idxs[1]);
        const auto& p2 = sys.position(idxs[2]);

        const auto r_ij         = sys.adjust_direction(p1, p0);
        const auto inv_len_r_ij = math::rlength(r_ij);
        const auto r_ij_reg     = r_ij * inv_len_r_ij;

        const auto r_kj         = sys.adjust_direction(p1, p2);
        const auto inv_len_r_kj = math::rlength(r_kj);
        const auto r_kj_reg     = r_kj * inv_len_r_kj;

        const auto dot_ijk   = math::dot_product(r_ij_reg, r_kj_reg);
        const auto cos_theta = math::clamp<real_type>(dot_ijk, -1, 1);

        // acos returns a value in [0, pi]
        const auto theta = std::acos(cos_theta);
        const auto coef  = -pot.derivative(theta);

        if(NeedEnergy)
        {
            energy += pot.potential(theta);
        }

        // sin(x) >= 0 if x is in [0, pi].
        const auto sin_theta    = std::sin(theta);
        const auto coef_inv_sin = coef / std::max(sin_theta, abs_tol);

        const auto Fi = (coef_inv_sin * inv_len_r_ij) *
                        (cos_theta * r_ij_reg - r_kj_reg);
        const auto Fk = (coef_inv_sin * inv_len_r_kj) *
                        (cos_theta * r_kj_reg - r_ij_reg);
        const auto Fj = -Fi - Fk;

        sys.force(idxs[0]) += Fi;
        sys.force(idxs[1]) += Fj;
        sys.force(idxs[2]) += Fk;

        if(NeedVirial)
        {
            sys.virial() += math::tensor_product(p1 + r_ij, Fi) +
                            math::tensor_product(p1,        Fj) +
                            math::tensor_product(p1 + r_kj, Fk);
        }
    }
    return energy;
}

template<typename traitsT, typename potentialT>
typename BondAngleInteraction<traitsT, potentialT>::real_type
BondAngleInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const auto& idxp : this->potentials_)
    {
        const auto& p0 = sys.position(idxp.first[0]);
        const auto& p1 = sys.position(idxp.first[1]);
        const auto& p2 = sys.position(idxp.first[2]);

        const auto rji = sys.adjust_direction(p1, p0); // pj -> pi
        const auto rjk = sys.adjust_direction(p1, p2); // pj -> pi

        const auto lsq_rji = math::length_sq(rji);
        const auto lsq_rjk = math::length_sq(rjk);
        const auto dot_ijk = math::dot_product(rji, rjk);

        const auto dot_ijk_reg = dot_ijk * math::rsqrt(lsq_rji * lsq_rjk);
        const auto cos_theta   = math::clamp<real_type>(dot_ijk_reg, -1, 1);
        const auto theta       = std::acos(cos_theta);

        E += idxp.second.potential(theta);
    }
    return E;
}

template<typename traitsT, typename potentialT>
void BondAngleInteraction<traitsT, potentialT>::write_topology(
        topology_type& topol) const
{
    if(this->kind_.empty() || this->kind_ == "none") {return;}

    for(const auto& idxp : this->potentials_)
    {
        const auto i = idxp.first[0];
        const auto j = idxp.first[1];
        const auto k = idxp.first[2];
        topol.add_connection(i, j, this->kind_);
        topol.add_connection(i, k, this->kind_);
        topol.add_connection(j, k, this->kind_);
    }
    return;
}

}// mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/HarmonicPotential.hpp>
#include <mjolnir/forcefield/local/GaussianPotential.hpp>
#include <mjolnir/forcefield/FLP/FlexibleLocalAnglePotential.hpp>

namespace mjolnir
{

// harmonic
extern template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, HarmonicPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, HarmonicPotential<float> >;
extern template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, HarmonicPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, HarmonicPotential<float> >;

// gaussian
extern template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, GaussianPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, GaussianPotential<float> >;
extern template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, GaussianPotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, GaussianPotential<float> >;

// FLP angle
extern template class BondAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, FlexibleLocalAnglePotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, FlexibleLocalAnglePotential<float> >;
extern template class BondAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<double>>;
extern template class BondAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, FlexibleLocalAnglePotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_BOND_ANGL_INTERACTION */
