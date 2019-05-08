#ifndef MJOLNIR_OMP_DIHEDRAL_ANGLE_INTERACTION_HPP
#define MJOLNIR_OMP_DIHEDRAL_ANGLE_INTERACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/interaction/local/DihedralAngleInteraction.hpp>

namespace mjolnir
{

// calculate energy and force of dihedral angle interaction.
// The implementation is based on the article written by H. Bekker,
// H.J.C. Berendsen and W.F. van Gunsteren (1995) J. Comput. Chem.
template<typename realT, template<typename, typename> class boundaryT,
         typename potentialT>
class DihedralAngleInteraction<OpenMPSimulatorTraits<realT, boundaryT>, potentialT>
    final : public LocalInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type          = OpenMPSimulatorTraits<realT, boundaryT>;
    using potential_type       = potentialT;
    using base_type            = LocalInteractionBase<traits_type>;
    using real_type            = typename base_type::real_type;
    using coordinate_type      = typename base_type::coordinate_type;
    using system_type          = typename base_type::system_type;
    using topology_type        = typename base_type::topology_type;
    using connection_kind_type = typename base_type::connection_kind_type;

    using indices_type         = std::array<std::size_t, 4>;
    using potential_index_pair = std::pair<indices_type, potentialT>;
    using container_type       = std::vector<potential_index_pair>;
    using iterator             = typename container_type::iterator;
    using const_iterator       = typename container_type::const_iterator;

  public:

    DihedralAngleInteraction(const connection_kind_type kind,
                             const container_type& pot)
        : kind_(kind), potentials(pot)
    {}
    DihedralAngleInteraction(const connection_kind_type kind,
                             container_type&& pot)
        : kind_(kind), potentials(std::move(pot))
    {}
    ~DihedralAngleInteraction() = default;

    void      calc_force (system_type& sys)       const noexcept override
    {
#pragma omp for nowait
        for(std::size_t i=0; i<this->potentials.size(); ++i)
        {
            const auto& idxp = this->potentials[i];
            const std::size_t idx0 = idxp.first[0];
            const std::size_t idx1 = idxp.first[1];
            const std::size_t idx2 = idxp.first[2];
            const std::size_t idx3 = idxp.first[3];

            const coordinate_type r_ij =
                sys.adjust_direction(sys.position(idx0) - sys.position(idx1));
            const coordinate_type r_kj =
                sys.adjust_direction(sys.position(idx2) - sys.position(idx1));
            const coordinate_type r_lk =
                sys.adjust_direction(sys.position(idx3) - sys.position(idx2));
            const coordinate_type r_kl = real_type(-1.0) * r_lk;

            const real_type r_kj_lensq  = math::length_sq(r_kj);
            const real_type r_kj_rlen   = math::rsqrt(r_kj_lensq);
            const real_type r_kj_rlensq = r_kj_rlen * r_kj_rlen;
            const real_type r_kj_len    = r_kj_rlen * r_kj_lensq;

            const coordinate_type m = math::cross_product(r_ij, r_kj);
            const coordinate_type n = math::cross_product(r_kj, r_kl);
            const real_type m_lensq = math::length_sq(m);
            const real_type n_lensq = math::length_sq(n);

            const coordinate_type R =
                r_ij - (math::dot_product(r_ij, r_kj) * r_kj_rlensq) * r_kj;
            const coordinate_type S =
                r_lk - (math::dot_product(r_lk, r_kj) * r_kj_rlensq) * r_kj;

            const real_type R_lensq = math::length_sq(R);
            const real_type S_lensq = math::length_sq(S);

            const real_type dot_RS  = math::dot_product(R, S) * math::rsqrt(R_lensq * S_lensq);
            const real_type cos_phi = math::clamp(dot_RS, real_type(-1.0), real_type(1.0));
            const real_type phi     =
                std::copysign(std::acos(cos_phi), math::dot_product(r_ij, n));

            // -dV / dphi
            const real_type coef = -(idxp.second.derivative(phi));

            const coordinate_type Fi = ( coef * r_kj_len / m_lensq) * m;
            const coordinate_type Fl = (-coef * r_kj_len / n_lensq) * n;

            const real_type coef_ijk = math::dot_product(r_ij, r_kj) * r_kj_rlensq;
            const real_type coef_jkl = math::dot_product(r_kl, r_kj) * r_kj_rlensq;

            sys.force_thread(omp_get_thread_num(), idx0) += Fi;
            sys.force_thread(omp_get_thread_num(), idx1) += (coef_ijk - real_type(1.0)) * Fi - coef_jkl * Fl;
            sys.force_thread(omp_get_thread_num(), idx2) += (coef_jkl - real_type(1.0)) * Fl - coef_ijk * Fi;
            sys.force_thread(omp_get_thread_num(), idx3) += Fl;
        }
        return;
    }

    real_type calc_energy(const system_type& sys) const noexcept override
    {
        real_type E = 0.0;
        for(const auto& idxp : this->potentials)
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

    void initialize(const system_type&) override
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        MJOLNIR_LOG_FUNCTION();
        MJOLNIR_LOG_INFO("With OpenMP: potential = ", potential_type::name(),
                         ", number of dihedrals = ", potentials.size());
        return;
    }

    void update(const system_type& sys) override
    {
        for(auto& item : potentials)
        {
            item.second.update(sys);
        }
    }

    // do nothing. this is used to reduce margin of neighbor list, and added
    // to this class for the consistency.
    void update_margin(const real_type, const system_type&) override {return;}

    std::string name() const override
    {return "DihedralAngle:"_s + potential_type::name();}

    void write_topology(topology_type& topol) const override
    {
        if(this->kind_.empty() || this->kind_ == "none") {return;}

        for(const auto& idxp : this->potentials)
        {
            const auto i = idxp.first[0];
            const auto j = idxp.first[1];
            const auto k = idxp.first[2];
            const auto l = idxp.first[3];
            topol.add_connection(i, j, this->kind_);
            topol.add_connection(j, k, this->kind_);
            topol.add_connection(k, l, this->kind_);
        }
        return;
    }

   private:
    connection_kind_type kind_;
    container_type potentials;
};

}// mjolnir
#endif /* MJOLNIR_DIHEDRAL_ANGLE_INTERACTION */
