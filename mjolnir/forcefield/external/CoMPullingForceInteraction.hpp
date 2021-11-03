#ifndef MJOLNIR_FORCEFIELD_EXTERNAL_COM_PULLING_FORCE_INTERACTION_HPP
#define MJOLNIR_FORCEFIELD_EXTERNAL_COM_PULLING_FORCE_INTERACTION_HPP
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/core/System.hpp>

namespace mjolnir
{

// pulling force interaction applies a constant force to the center of mass (com)
// of the specified atoms without any particular reason.
template<typename traitsT>
class CoMPullingForceInteraction final
    : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using base_type       = ExternalForceInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

    // parameter = pairof {{particle indices}, force}
    using parameter_type = std::pair<std::vector<std::size_t>, coordinate_type>;
    using container_type = std::vector<parameter_type>;

  public:

    CoMPullingForceInteraction(std::vector<parameter_type>&& parameters)
        : parameters_(std::move(parameters))
    {}
    ~CoMPullingForceInteraction() override {}

    void calc_force (system_type& sys) const noexcept override
    {
        assert(parameters_.size() == total_masses_.size());
        for(std::size_t idx=0; idx<parameters_.size(); ++idx)
        {
            const auto& para          = parameters_[idx];
            const auto total_mass     = total_masses_[idx];
            const auto inv_total_mass = real_type(1) / total_mass;
            const auto force          = para.second * inv_total_mass;
            for(const auto i : para.first)
            {
                sys.force(i) += sys.mass(i) * force;
            }
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

    void initialize(const system_type& sys) override
    {
        this->update(sys);
        return;
    }
    void update(const system_type&) override
    {
        total_masses_.clear();
        total_masses_.reserve(parameters_.size());
        for(const auto& para : parameters_)
        {
            real_type total_mass(0);
            for(const auto i : para.first)
            {
                total_mass += sys.mass(i);
            }
            total_masses_.push_back(total_mass);
        }
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

    std::string name() const override {return std::string("CoMPullingForce:");}

    base_type* clone() const override
    {
        return new CoMPullingForceInteraction(*this);
    }

    container_type const& parameters() const noexcept {return parameters_;}

  private:

    std::vector<real_type> total_masses_;
    container_type parameters_;
};

} // mjolnir
#endif // MJOLNIR_FORCEFIELD_EXTERNAL_PULLING_FORCE_INTERACTION_HPP
