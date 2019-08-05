#ifndef MJOLNIR_POTENTIAL_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>

namespace mjolnir
{

/*! @brief Lennard-Jones form wall potential                               *
 * designed for external force field. For particle-particle L-J potential, *
 * see mjolnir/potential/global/LennardJonesPotential.hpp.                 */
template<typename realT>
class LennardJonesWallPotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::pair<real_type, real_type>;

  public:

    LennardJonesWallPotential(
        const real_type cutoff_ratio, std::vector<parameter_type> params)
        : cutoff_ratio_(cutoff_ratio),
          coef_at_cutoff_(std::pow(1 / cutoff_ratio, 12) -
                          std::pow(1 / cutoff_ratio, 6)),
          params_(std::move(params))
    {}
    LennardJonesWallPotential(std::vector<parameter_type> params)
        : LennardJonesWallPotential(2.5, std::move(params))
    {}
    ~LennardJonesWallPotential(){}

    real_type potential(const std::size_t i, const real_type r) const noexcept
    {
        const real_type sigma = params_[i].first;
        if(sigma * this->cutoff_ratio_ < r){return real_type(0.0);}
        const real_type epsilon = params_[i].second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return real_type(4.0) * epsilon * (r12s12 - r6s6 - this->coef_at_cutoff_);
    }

    real_type derivative(const std::size_t i, const real_type r) const noexcept
    {
        const real_type sigma = params_[i].first;
        if(sigma * this->cutoff_ratio_ < r){return real_type(0.0);}
        const real_type epsilon = params_[i].second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return real_type(24.0) * epsilon * (r6s6 - real_type(2.0) * r12s12) / r;
    }

    real_type max_cutoff_length() const noexcept
    {
        const real_type max_sigma = std::max_element(
            this->params_.cbegin(), this->params_.cend(),
            [](const parameter_type& lhs, const parameter_type& rhs) noexcept {
                return lhs.first < rhs.first;
            })->first;
        return max_sigma * this->cutoff_ratio_;
    }

    std::vector<std::size_t> participants() const
    {
        std::vector<std::size_t> retval;
        retval.reserve(this->params_.size());
        for(std::size_t i=0; i<this->params_.size(); ++i)
        {
            if(this->params_[i].second != real_type(0.0))
            {
                retval.push_back(i);
            }
        }
        return retval;
    }

    // nothing to do when system parameters change.
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    const char* name() const noexcept {return "LennardJonesWall";}

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return params_;}
    std::vector<parameter_type> const& parameters() const noexcept {return params_;}

  private:

    real_type cutoff_ratio_, coef_at_cutoff_;
    std::vector<parameter_type> params_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class LennardJonesWallPotential<double>;
extern template class LennardJonesWallPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif // MJOLNIR_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL
