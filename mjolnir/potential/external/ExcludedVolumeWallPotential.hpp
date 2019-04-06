#ifndef MJOLNIR_EXTERNAL_EXCLUDED_VOLUME_WALL_POTENTIAL
#define MJOLNIR_EXTERNAL_EXCLUDED_VOLUME_WALL_POTENTIAL
#include <mjolnir/core/System.hpp>
#include <mjolnir/math/math.hpp>

namespace mjolnir
{

/*! @brief excluded volume potential                                       *
 * designed for external force field. For particle-particle EXV potential, *
 * see mjolnir/potential/global/ExcludedVolumePotential.hpp.               */
template<typename realT>
class ExcludedVolumeWallPotential
{
  public:
    using real_type = realT;

    // rc = 2.0 * sigma
    constexpr static real_type cutoff_ratio = 2.0;
    // to make the potential curve continuous at the cutoff point
    constexpr static real_type coef_at_cutoff =
        compiletime::pow(1.0 / cutoff_ratio, 12);

  public:

    ExcludedVolumeWallPotential(
        const real_type epsilon, std::vector<real_type> params)
        : epsilon_(epsilon), radii_(std::move(params))
    {}
    ~ExcludedVolumeWallPotential(){}

    real_type potential(const std::size_t i, const real_type r) const noexcept
    {
        const real_type d = this->radii_[i];
        if(d * cutoff_ratio < r){return real_type(0.0);}

        const real_type d_r  = d / r;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return this->epsilon_ * (dr12 - coef_at_cutoff);
    }

    real_type derivative(const std::size_t i, const real_type r) const noexcept
    {
        const real_type d = this->radii_[i];
        if(d * cutoff_ratio < r){return real_type(0.0);}

        const real_type rinv = real_type(1.0) / r;
        const real_type d_r  = d * rinv;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return real_type(-12.0) * this->epsilon_ * dr12 * rinv;
    }

    real_type max_cutoff_length() const noexcept
    {
        const real_type max_sigma =
            *(std::max_element(radii_.cbegin(), radii_.cend()));
        return max_sigma * cutoff_ratio;
    }

    // TODO more sophisticated way to treat this
    std::vector<std::size_t> participants() const
    {
        std::vector<std::size_t> retval;
        retval.reserve(this->radii_.size());
        for(std::size_t i=0; i<this->radii_.size(); ++i)
        {
            if(this->radii_[i] != real_type(0.0))
            {
                retval.push_back(i);
            }
        }
        return retval;
    }

    // nothing to do when system parameters change.
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    const char* name() const noexcept {return "ExcludedVolumeWall";}

    real_type epsilon() const noexcept {return epsilon_;}

    // access to the parameters...
    std::vector<real_type>&       parameters()       noexcept {return radii_;}
    std::vector<real_type> const& parameters() const noexcept {return radii_;}

  private:

    real_type epsilon_;
    std::vector<real_type> radii_;
};

template<typename realT>
constexpr typename ExcludedVolumeWallPotential<realT>::real_type
ExcludedVolumeWallPotential<realT>::cutoff_ratio;
template<typename realT>
constexpr typename ExcludedVolumeWallPotential<realT>::real_type
ExcludedVolumeWallPotential<realT>::coef_at_cutoff;

} // mjolnir
#endif // MJOLNIR_EXTERNAL_EXCLUDED_VOLUME_WALL_POTENTIAL
