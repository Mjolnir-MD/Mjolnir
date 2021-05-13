#ifndef MJOLNIR_FORCEFIELD_EXTERNAL_PULLING_FORCE_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_EXTERNAL_PULLING_FORCE_INTERACTION_HPP
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

// pulling force interaction applies a constant force to specified atoms
// without any particular reason.
template<typename traitsT>
class PullingForceInteraction final
    : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = ExternalForceInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

    using parameter_type = std::pair<std::size_t, coordinate_type>;// {idx, force}
    using container_type = std::vector<parameter_type>;

  public:

    PullingForceInteraction(container_type&& parameters)
        : parameters_(std::move(parameters))
    {}
    ~PullingForceInteraction() override {}

    // calculate force, update spatial partition (reduce margin) inside.
    void      calc_force (system_type& sys)           const noexcept override
    {
        for(const auto& param : parameters_)
        {
            sys.force(param.first) += param.second;
        }
        return;
    }
    real_type calc_energy(system_type const&)         const noexcept override
    {
        return real_type(0);
    }
    real_type calc_force_and_energy(system_type& sys) const noexcept override
    {
        this->calc_force(sys);
        return this->calc_energy(sys);
    }

    void initialize(const system_type&) override
    {
        return;
    }
    void update(const system_type&) override
    {
        return;
    }

    void reduce_margin(const real_type, const system_type&) override
    {
        return;
    }
    void scale_margin(const real_type, const system_type&) override
    {
        return;
    }

    std::string name() const override {return "PullingForce:"_s;}

    base_type* clone() const override
    {
        return new PullingForceInteraction(container_type(parameters_));
    }

  private:

    container_type parameters_;
};

} // mjolnir
#endif // MJOLNIR_FORCEFIELD_EXTERNAL_PULLING_FORCE_INTERACTION_HPP
