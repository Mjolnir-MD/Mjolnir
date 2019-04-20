#ifndef MJOLNIR_POTENTIAL_GLOBAL_UNIFORM_LENNARD_JONES_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_GLOBAL_UNIFORM_LENNARD_JONES_POTENTIAL_HPP
#include <mjolnir/potential/global/IgnoreMolecule.hpp>
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/empty.hpp>
#include <vector>
#include <algorithm>
#include <cmath>

namespace mjolnir
{

// Well-known Lennard-Jones interaction with uniform parameters.
// This class contains a sigma and an epsilon that are the same among all the
// particles.
template<typename realT>
class UniformLennardJonesPotential
{
  public:
    using real_type = realT;
    using parameter_type = empty_t;
    using container_type = empty_t;

    // topology stuff
    using topology_type        = Topology;
    using molecule_id_type     = typename topology_type::molecule_id_type;
    using connection_kind_type = typename topology_type::connection_kind_type;
    using ignore_molecule_type = IgnoreMolecule<molecule_id_type>;

    // rc = 2.5 * sigma
    constexpr static real_type cutoff_ratio = 2.5;
    // to make the potential curve continuous at the cutoff point
    constexpr static real_type coef_at_cutoff =
        compiletime::pow<real_type>(1 / cutoff_ratio, 12u) -
        compiletime::pow<real_type>(1 / cutoff_ratio,  6u);

  public:

    UniformLennardJonesPotential(const real_type sgm, const real_type eps,
        const std::map<connection_kind_type, std::size_t>& exclusions,
        ignore_molecule_type ignore_molecule)
        : sigma_(sgm), epsilon_(eps), r_cut_(sgm * cutoff_ratio),
          ignore_molecule_(std::move(ignore_molecule)),
          ignore_within_(exclusions.begin(), exclusions.end())
    {}
    ~UniformLennardJonesPotential() = default;

    parameter_type prepare_params(std::size_t, std::size_t) const noexcept
    {
        return parameter_type{}; // no pre-calculated parameter
    }

    // forwarding functions for clarity...
    real_type potential(const std::size_t i, const std::size_t j,
                        const real_type r) const noexcept
    {
        return this->potential(r, this->prepare_params(i, j));
    }
    real_type derivative(const std::size_t i, const std::size_t j,
                         const real_type r) const noexcept
    {
        return this->derivative(r, this->prepare_params(i, j));
    }

    real_type
    potential(const real_type r, const parameter_type&) const noexcept
    {
        if(r_cut_ < r){return 0;}

        const real_type r1s1   = sigma_ / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 4 * epsilon_ * (r12s12 - r6s6 - coef_at_cutoff);
    }

    real_type
    derivative(const real_type r, const parameter_type&) const noexcept
    {
        if(r_cut_ < r){return 0;}

        const real_type rinv   = 1 / r;
        const real_type r1s1   = sigma_ * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return 24 * epsilon_ * (r6s6 - 2 * r12s12) * rinv;
    }

    real_type max_cutoff_length() const noexcept
    {
        return sigma_ * cutoff_ratio;
    }

    // nothing to do when system parameters change.
    template<typename traitsT>
    void update(const System<traitsT>&) const noexcept {return;}

    // e.g. `3` means ignore particles connected within 3 "bond"s
    std::vector<std::pair<connection_kind_type, std::size_t>>
    ignore_within() const {return ignore_within_;}

    bool is_ignored_molecule(
            const molecule_id_type& i, const molecule_id_type& j) const noexcept
    {
        return ignore_molecule_.is_ignored(i, j);
    }

    static const char* name() noexcept {return "LennardJones";}

    // access to the parameters...
    real_type& sigma()         noexcept {return sigma_;}
    real_type  sigma()   const noexcept {return sigma_;}
    real_type& epsilon()       noexcept {return epsilon_;}
    real_type  epsilon() const noexcept {return epsilon_;}

  private:

    real_type sigma_, epsilon_, r_cut_;

    ignore_molecule_type ignore_molecule_;
    std::vector<std::pair<connection_kind_type, std::size_t>> ignore_within_;
};

template<typename traitsT>
constexpr typename UniformLennardJonesPotential<traitsT>::real_type
UniformLennardJonesPotential<traitsT>::cutoff_ratio;

} // mjolnir
#endif /* MJOLNIR_LENNARD_JONES_POTENTIAL */
