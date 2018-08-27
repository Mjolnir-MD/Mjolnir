#ifndef MJOLNIR_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL
#define MJOLNIR_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL
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

    // rc = 2.5 * sigma
    constexpr static real_type cutoff_ratio = 2.5;
    // to make the potential curve continuous at the cutoff point
    constexpr static real_type coef_at_cutoff =
        compiletime::pow(1.0 / cutoff_ratio, 12u) -
        compiletime::pow(1.0 / cutoff_ratio,  6u);

  public:

    LennardJonesWallPotential(std::vector<parameter_type> params)
        : params_(std::move(params))
    {}
    ~LennardJonesWallPotential(){}

    real_type potential(const std::size_t i, const real_type r) const noexcept
    {
        const real_type sigma = params_[i].first;
        if(sigma * cutoff_ratio < r){return real_type(0.0);}
        const real_type epsilon = params_[i].second;

        const real_type r1s1   = sigma / r;
        const real_type r3s3   = r1s1 * r1s1 * r1s1;
        const real_type r6s6   = r3s3 * r3s3;
        const real_type r12s12 = r6s6 * r6s6;
        return real_type(4.0) * epsilon * (r12s12 - r6s6 - coef_at_cutoff);
    }

    real_type derivative(const std::size_t i, const real_type r) const noexcept
    {
        const real_type sigma = params_[i].first;
        if(sigma * cutoff_ratio < r){return real_type(0.0);}
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
        return max_sigma * cutoff_ratio;
    }

    // TODO
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
    void update(const System<T>& sys) const noexcept {return;}

    const char* name() const noexcept {return "LennardJonesWall";}

    // access to the parameters...
    std::vector<parameter_type>&       params()       noexcept {return params_;}
    std::vector<parameter_type> const& params() const noexcept {return params_;}

  private:

    std::vector<parameter_type> params_;
};

template<typename realT>
constexpr typename LennardJonesWallPotential<realT>::real_type
LennardJonesWallPotential<realT>::cutoff_ratio;
template<typename realT>
constexpr typename LennardJonesWallPotential<realT>::real_type
LennardJonesWallPotential<realT>::coef_at_cutoff;

} // mjolnir
#endif // MJOLNIR_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL
