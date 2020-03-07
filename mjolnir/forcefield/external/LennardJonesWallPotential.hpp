#ifndef MJOLNIR_POTENTIAL_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>

namespace mjolnir
{

/*! @brief Lennard-Jones form wall potential                               *
 * designed for external force field. For particle-particle L-J potential, *
 * see mjolnir/forcefield/global/LennardJonesPotential.hpp.                 */
template<typename realT>
class LennardJonesWallPotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::pair<real_type, real_type>;
    using container_type = std::vector<parameter_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2.5);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return parameter_type{real_type(0), real_type(0)};
    }

  public:

    LennardJonesWallPotential(const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, parameter_type>>& params)
        : cutoff_ratio_(cutoff_ratio), coef_at_cutoff_(
            std::pow(1 / cutoff_ratio, 12) - std::pow(1 / cutoff_ratio, 6))
    {
        this->parameters_.resize(params.size());
        this->participants_.reserve(params.size());
        for(const auto& idxp: params)
        {
            const auto idx = idxp.first;
            this->participants_.push_back(idx);
            if(idx >= this->parameters_.size())
            {
                this->parameters_.resize(idx+1, default_parameter());
            }
            this->parameters_.at(idx) = idxp.second;
        }
    }
    ~LennardJonesWallPotential(){}

    real_type potential(const std::size_t i, const real_type r) const noexcept
    {
        const real_type sigma = parameters_[i].first;
        if(sigma * this->cutoff_ratio_ < r){return real_type(0.0);}
        const auto eps    = parameters_[i].second;
        const auto r1s1   = sigma / r;
        const auto r3s3   = r1s1 * r1s1 * r1s1;
        const auto r6s6   = r3s3 * r3s3;
        const auto r12s12 = r6s6 * r6s6;
        return real_type(4.0) * eps * (r12s12 - r6s6 - this->coef_at_cutoff_);
    }
    real_type derivative(const std::size_t i, const real_type r) const noexcept
    {
        const real_type sigma = parameters_[i].first;
        if(sigma * this->cutoff_ratio_ < r){return real_type(0.0);}

        const real_type eps    = parameters_[i].second;
        const real_type rinv   = real_type(1) / r;
        const real_type r1s1   = sigma * rinv;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return real_type(24.0) * eps * (r6s6 - real_type(2.0) * r12s12) * rinv;
    }
    real_type max_cutoff_length() const noexcept
    {
        const real_type max_sigma = std::max_element(
            this->parameters_.cbegin(), this->parameters_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return lhs.first < rhs.first;
            })->first;
        return max_sigma * this->cutoff_ratio_;
    }

    std::vector<std::size_t> const& participants() const noexcept
    {
        return participants_;
    }

    // nothing to do when system parameters change.
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    const char* name() const noexcept {return "LennardJonesWall";}

    // access to the parameters...
    container_type&       parameters()       noexcept {return parameters_;}
    container_type const& parameters() const noexcept {return parameters_;}

    real_type cutoff() const noexcept {return this->cutoff_ratio_;}

  private:

    real_type cutoff_ratio_, coef_at_cutoff_;
    container_type           parameters_;
    std::vector<std::size_t> participants_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class LennardJonesWallPotential<double>;
extern template class LennardJonesWallPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif // MJOLNIR_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL
