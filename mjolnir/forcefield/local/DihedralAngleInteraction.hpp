#ifndef MJOLNIR_INTERACTION_DIHEDRAL_ANGLE_INTERACTION_HPP
#define MJOLNIR_INTERACTION_DIHEDRAL_ANGLE_INTERACTION_HPP
#include <mjolnir/core/LocalInteractionBase.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/string.hpp>
#include <mjolnir/util/logger.hpp>
#include <cmath>

namespace mjolnir
{

// calculate energy and force of dihedral angle interaction.
// The implementation is based on the article written by H. Bekker,
// H.J.C. Berendsen and W.F. van Gunsteren (1995) J. Comput. Chem.
template<typename traitsT, typename potentialT>
class DihedralAngleInteraction final : public LocalInteractionBase<traitsT>
{
  public:
    using traits_type           = traitsT;
    using potential_type        = potentialT;
    using base_type             = LocalInteractionBase<traits_type>;
    using real_type             = typename base_type::real_type;
    using coordinate_type       = typename base_type::coordinate_type;
    using system_type           = typename base_type::system_type;
    using topology_type         = typename base_type::topology_type;
    using connection_kind_type  = typename base_type::connection_kind_type;

    using indices_type          = std::array<std::size_t, 4>;
    using potential_index_pair  = std::pair<indices_type, potentialT>;
    using container_type        = std::vector<potential_index_pair>;
    using iterator              = typename container_type::iterator;
    using const_iterator        = typename container_type::const_iterator;

  public:

    DihedralAngleInteraction(const connection_kind_type kind,
                             const container_type& pot)
        : kind_(kind), potentials_(pot)
    {}
    DihedralAngleInteraction(const connection_kind_type kind,
                             container_type&& pot)
        : kind_(kind), potentials_(std::move(pot))
    {}
    ~DihedralAngleInteraction() {}

    void      calc_force (system_type&)       const noexcept override;
    real_type calc_energy(const system_type&) const noexcept override;

    void initialize(const system_type& sys) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("potential = ", potential_type::name(),
                         ", number of dihedrals = ", potentials_.size());
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
    {return "DihedralAngle:"_s + potential_type::name();}

    void write_topology(topology_type&) const override;

    container_type const& potentials() const noexcept {return potentials_;}
    container_type&       potentials()       noexcept {return potentials_;}

    base_type* clone() const override
    {
        return new DihedralAngleInteraction(kind_, container_type(potentials_));
    }

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


template<typename traitsT, typename pT>
void
DihedralAngleInteraction<traitsT, pT>::calc_force(system_type& sys) const noexcept
{
    for(const auto& idxp : this->potentials_)
    {
        const std::size_t idx0 = idxp.first[0];
        const std::size_t idx1 = idxp.first[1];
        const std::size_t idx2 = idxp.first[2];
        const std::size_t idx3 = idxp.first[3];

        const auto& r_i = sys.position(idx0);
        const auto& r_j = sys.position(idx1);
        const auto& r_k = sys.position(idx2);
        const auto& r_l = sys.position(idx3);

        const coordinate_type r_ij = sys.adjust_direction(r_i - r_j);
        const coordinate_type r_kj = sys.adjust_direction(r_k - r_j);
        const coordinate_type r_lk = sys.adjust_direction(r_l - r_k);
        const coordinate_type r_kl = real_type(-1.0) * r_lk;

        const real_type r_kj_lensq  = math::length_sq(r_kj);
        const real_type r_kj_rlen   = math::rsqrt(r_kj_lensq);
        const real_type r_kj_rlensq = r_kj_rlen * r_kj_rlen;
        const real_type r_kj_len    = r_kj_rlen * r_kj_lensq;

        const coordinate_type m = math::cross_product(r_ij, r_kj);
        const coordinate_type n = math::cross_product(r_kj, r_kl);
        const real_type m_lensq = math::length_sq(m);
        const real_type n_lensq = math::length_sq(n);

        const real_type dot_mn  = math::dot_product(m, n) *
                                  math::rsqrt(m_lensq * n_lensq);
        const real_type cos_phi = math::clamp<real_type>(dot_mn, -1, 1);
        const real_type phi     =
            std::copysign(std::acos(cos_phi), math::dot_product(r_ij, n));

        // -dV / dphi
        const real_type coef = -(idxp.second.derivative(phi));

        const coordinate_type Fi = ( coef * r_kj_len / m_lensq) * m;
        const coordinate_type Fl = (-coef * r_kj_len / n_lensq) * n;

        const real_type coef_ijk = math::dot_product(r_ij, r_kj) * r_kj_rlensq;
        const real_type coef_jkl = math::dot_product(r_kl, r_kj) * r_kj_rlensq;

        sys.force(idx0) += Fi;
        sys.force(idx1) += (coef_ijk - real_type(1.0)) * Fi - coef_jkl * Fl;
        sys.force(idx2) += (coef_jkl - real_type(1.0)) * Fl - coef_ijk * Fi;
        sys.force(idx3) += Fl;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename DihedralAngleInteraction<traitsT, potentialT>::real_type
DihedralAngleInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const auto& idxp : this->potentials_)
    {
        const coordinate_type r_ij = sys.adjust_direction(
                sys.position(idxp.first[0]) - sys.position(idxp.first[1]));
        const coordinate_type r_kj = sys.adjust_direction(
                sys.position(idxp.first[2]) - sys.position(idxp.first[1]));
        const coordinate_type r_lk = sys.adjust_direction(
                sys.position(idxp.first[3]) - sys.position(idxp.first[2]));
        const real_type r_kj_lensq_inv = real_type(1.0) / math::length_sq(r_kj);

        const coordinate_type n = math::cross_product(r_kj, real_type(-1.0) * r_lk);

        const coordinate_type R = r_ij -
                              (math::dot_product(r_ij, r_kj) * r_kj_lensq_inv) * r_kj;
        const coordinate_type S = r_lk -
                              (math::dot_product(r_lk, r_kj) * r_kj_lensq_inv) * r_kj;
        const real_type R_lensq = math::length_sq(R);
        const real_type S_lensq = math::length_sq(S);

        const real_type dot_RS  = math::dot_product(R, S) * math::rsqrt(R_lensq * S_lensq);
        const real_type cos_phi = math::clamp(dot_RS, real_type(-1.0), real_type(1.0));
        const real_type phi     =
            std::copysign(std::acos(cos_phi), math::dot_product(r_ij, n));

        E += idxp.second.potential(phi);
    }
    return E;
}

template<typename traitsT, typename potentialT>
void DihedralAngleInteraction<traitsT, potentialT>::write_topology(
        topology_type& topol) const
{
    if(this->kind_.empty() || this->kind_ == "none") {return;}

    for(const auto& idxp : this->potentials_)
    {
        const auto i = idxp.first[0];
        const auto j = idxp.first[1];
        const auto k = idxp.first[2];
        const auto l = idxp.first[3];
        topol.add_connection(i, j, this->kind_);
        topol.add_connection(i, k, this->kind_);
        topol.add_connection(i, l, this->kind_);
        topol.add_connection(j, k, this->kind_);
        topol.add_connection(j, l, this->kind_);
        topol.add_connection(k, l, this->kind_);
    }
    return;
}

}// mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/forcefield/local/ClementiDihedralPotential.hpp>
#include <mjolnir/forcefield/local/PeriodicGaussianPotential.hpp>
#include <mjolnir/forcefield/local/CosinePotential.hpp>
#include <mjolnir/forcefield/FLP/FlexibleLocalDihedralPotential.hpp>

namespace mjolnir
{

// ClementiDihedral
extern template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, ClementiDihedralPotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, ClementiDihedralPotential<float> >;
extern template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ClementiDihedralPotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ClementiDihedralPotential<float> >;

// gaussian
extern template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, PeriodicGaussianPotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, PeriodicGaussianPotential<float> >;
extern template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, PeriodicGaussianPotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, PeriodicGaussianPotential<float> >;

// FLP dihedral
extern template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, FlexibleLocalDihedralPotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, FlexibleLocalDihedralPotential<float> >;
extern template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, FlexibleLocalDihedralPotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, FlexibleLocalDihedralPotential<float> >;

// cosine
extern template class DihedralAngleInteraction<SimulatorTraits<double, UnlimitedBoundary>, CosinePotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  UnlimitedBoundary>, CosinePotential<float> >;
extern template class DihedralAngleInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, CosinePotential<double>>;
extern template class DihedralAngleInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, CosinePotential<float> >;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif /* MJOLNIR_DIHEDRAL_ANGLE_INTERACTION */
