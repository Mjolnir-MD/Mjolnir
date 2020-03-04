#ifndef MJOLNIR_POTENTIAL_EXTERNAL_EXCLUDED_VOLUME_WALL_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_EXTERNAL_EXCLUDED_VOLUME_WALL_POTENTIAL_HPP
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
    using parameter_type = real_type;
    using container_type = std::vector<parameter_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return real_type(2);
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return real_type(0);
    }

  public:

    ExcludedVolumeWallPotential(
        const real_type epsilon, const real_type cutoff_ratio,
        const std::vector<std::pair<std::size_t, real_type>>& params)
        : epsilon_(epsilon), cutoff_ratio_(cutoff_ratio),
          coef_at_cutoff_(std::pow(1 / cutoff_ratio, 12))
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
    ~ExcludedVolumeWallPotential(){}

    real_type potential(const std::size_t i, const real_type r) const noexcept
    {
        const real_type d = this->parameters_[i];
        if(d * this->cutoff_ratio_ < r){return real_type(0.0);}

        const real_type d_r  = d / r;
        const real_type dr3  = d_r * d_r * d_r;
        const real_type dr6  = dr3 * dr3;
        const real_type dr12 = dr6 * dr6;
        return this->epsilon_ * (dr12 - this->coef_at_cutoff_);
    }

    real_type derivative(const std::size_t i, const real_type r) const noexcept
    {
        const real_type d = this->parameters_[i];
        if(d * this->cutoff_ratio_ < r){return real_type(0.0);}

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
            *(std::max_element(parameters_.cbegin(), parameters_.cend()));
        return max_sigma * this->cutoff_ratio_;
    }

    std::vector<std::size_t> const& participants() const noexcept
    {
        return participants_;
    }

    // nothing to do when system parameters change.
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    const char* name() const noexcept {return "ExcludedVolumeWall";}

    real_type epsilon() const noexcept {return epsilon_;}

    // access to the parameters...
    container_type&       parameters()       noexcept {return parameters_;}
    container_type const& parameters() const noexcept {return parameters_;}

    real_type cutoff() const noexcept {return this->cutoff_ratio_;}

  private:

    real_type epsilon_;
    real_type cutoff_ratio_, coef_at_cutoff_;
    container_type           parameters_;
    std::vector<std::size_t> participants_;
};

#ifdef MJOLNIR_SEPARATE_BUILD
extern template class ExcludedVolumeWallPotential<double>;
extern template class ExcludedVolumeWallPotential<float>;
#endif// MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif // MJOLNIR_EXTERNAL_EXCLUDED_VOLUME_WALL_POTENTIAL
