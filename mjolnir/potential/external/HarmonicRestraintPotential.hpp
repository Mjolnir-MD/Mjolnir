#ifndef MJOLNIR_POTENTIAL_EXTERNAL_HARMONIC_RESTRAINT_POTENTIAL_HPP
#define MJOLNIR_POTENTIAL_EXTERNAL_HARMONIC_RESTRAINT_POTENTIAL_HPP
#include <mjolnir/core/System.hpp>
#include <utility>
#include <vector>
#include <limits>

namespace mjolnir
{

// apply a harmonic restraint (k(x-x0)^2.) per one particle.
template<typename realT>
class HarmonicRestraintPotential
{
  public:
    using real_type      = realT;
    using parameter_type = std::pair<real_type, real_type>; // {k, x0}

  public:

    HarmonicRestraintPotential(std::vector<parameter_type> params)
        : params_(std::move(params))
    {}
    ~HarmonicRestraintPotential() = default;

    real_type potential(const std::size_t i, const real_type r) const noexcept
    {
        const real_type k  = this->params_[i].first;
        const real_type r0 = this->params_[i].second;
        const real_type dr = r - r0;

        return k * dr * dr;
    }

    real_type derivative(const std::size_t i, const real_type r) const noexcept
    {
        const real_type k  = this->params_[i].first;
        const real_type r0 = this->params_[i].second;
        const real_type dr = r - r0;

        return 2 * k * dr;
    }

    // no cutoff length. it has a finite value wherever particle is.
    real_type max_cutoff_length() const noexcept
    {
        return std::numeric_limits<real_type>::infinity();
    }

    std::vector<std::size_t> participants() const
    {
        std::vector<std::size_t> retval;
        retval.reserve(this->params_.size());
        for(std::size_t i=0; i<this->params_.size(); ++i)
        {
            if(this->params_[i].first != real_type(0.0))
            {
                retval.push_back(i);
            }
        }
        return retval;
    }

    // nothing to do when system parameters change.
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    const char* name() const noexcept {return "HarmonicRestraint";}

    // access to the parameters...
    std::vector<parameter_type>&       parameters()       noexcept {return params_;}
    std::vector<parameter_type> const& parameters() const noexcept {return params_;}

  private:

    std::vector<parameter_type> params_;
};

} // mjolnir
#endif // MJOLNIR_EXTERNAL_LENNARD_JONES_WALL_POTENTIAL
