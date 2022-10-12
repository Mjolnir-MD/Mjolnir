#ifndef MJOLNIR_POTENTIAL_EXTERNAL_HARMONIC_GROOVE_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_EXTERNAL_HARMONIC_GROOVE_POTENTIAL_HPP
#include <mjolnir/core/System.hpp>
#include <limits>

namespace mjolnir
{

/*! @brief harmonic potential                                                     *
 * designed for external force field. Force particle-particle Harmonic potential, *
 * see mjolnir/forcefield/local/HarmonicPotential.hpp                             */
template<typename realT>
class HarmonicGroovePotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::pair<real_type, real_type>;
    using container_type = std::vector<parameter_type>;

    static constexpr real_type default_cutoff() noexcept
    {
        return std::numeric_limits<real_type>::max();
    }
    static constexpr parameter_type default_parameter() noexcept
    {
        return std::make_pair(real_type(0), real_type(0));
    }

  public:

    HarmonicGroovePotential(
        const std::vector<std::pair<std::size_t, parameter_type>>& params)
    {
        this->parameters_.resize(params.size());
        this->participants_.reserve(params.size());
        for(const auto& idxp : params)
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
    ~HarmonicGroovePotential(){}

    real_type potential(const std::size_t i, const real_type v) const noexcept
    {
        const parameter_type k_v0 = this->parameters_[i];
        const real_type dv   = v - k_v0.second;
        return k_v0.first * dv * dv;
    }

    real_type derivative(const std::size_t i, const real_type v) const noexcept
    {
        const parameter_type k_v0 = this->parameters_[i];
        return 2 * k_v0.first * (v - k_v0.second);
    }

    real_type max_cutoff_length() const noexcept
    {
        return this->default_cutoff();
    }

    std::vector<std::size_t> const& participants() const noexcept
    {
        return participants_;
    }

    // nothing to do when system parameters change.
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    const char* name() const noexcept {return "HarmonicGroove";}

    container_type&       parameters()       noexcept {return parameters_;}
    container_type const& parameters() const noexcept {return parameters_;}

  private:

    container_type           parameters_;
    std::vector<std::size_t> participants_;
};

#ifndef MJOLNIR_SEPARATE_BUILD
extern template class HarmonicGroovePotential<double>;
extern template class HarmonicGroovePotential<float>;
#endif // MJOLNIR_SEPARATE_BUILD

} // mjolnir
#endif // MJOLNIR_POTENTIAL_EXTERNAL_HARMONIC_GROOVE_POTENTIAL_HPP
