#ifndef MJOLNIR_OMP_EXTERNAL_POSITION_RESTRAINT_INTERACTION_HPP
#define MJOLNIR_OMP_EXTERNAL_POSITION_RESTRAINT_INTERACTION_HPP
#include <mjolnir/omp/OpenMPSimulatorTraits.hpp>
#include <mjolnir/interaction/external/PositionRestraintInteraction.hpp>

namespace mjolnir
{

//! @brief Interaction between particles and points.
template<typename realT, template<typename, typename> class boundaryT,
         typename potentialT>
class PositionRestraintInteraction<
    OpenMPSimulatorTraits<realT, boundaryT>, potentialT
    > final: public ExternalForceInteractionBase<OpenMPSimulatorTraits<realT, boundaryT>>
{
  public:
    using traits_type     = OpenMPSimulatorTraits<realT, boundaryT>;
    using potential_type  = potentialT;
    using base_type       = ExternalForceInteractionBase<traits_type>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

    // list of {pair of a particle-id and a position at which
    //          the particle would be restrained}
    using shape_type = std::vector<std::pair<std::size_t, coordinate_type>>;

  public:

    PositionRestraintInteraction(
        shape_type shape, potential_type pot)
        : shape_(std::move(shape)), potential_(std::move(pot))
    {}
    ~PositionRestraintInteraction() override = default;

    void calc_force (system_type& sys) const noexcept override
    {
#pragma omp for nowait
        for(std::size_t i=0; i<this->shape_.size(); ++i)
        {
            const auto& pidpos = this->shape_[i];

            const auto  pid = pidpos.first;  // particle index
            const auto& pos = pidpos.second; // restrained position

            const auto& r_i = sys.position(pid);
            const auto  dr  = sys.adjust_direction(pos - r_i);

            const auto rlen = math::rlength(dr); // 1 / |dr|
            const auto dist = math::length_sq(dr) * rlen;
            const auto dV   = this->potential_.derivative(pid, dist);
            if(dV == 0.0){continue;}

            sys.force_thread(omp_get_thread_num(), pid) += (dV * rlen) * dr;
        }
        return ;
    }

    real_type calc_energy(system_type const& sys) const noexcept override
    {
        real_type E = 0.0;
#pragma omp parallel for reduction(+:E)
        for(std::size_t i=0; i<this->shape_.size(); ++i)
        {
            const auto& pidpos = this->shape_[i];
            const auto  pid = pidpos.first;  // particle index
            const auto& pos = pidpos.second; // restrained position

            const auto& r_i = sys.position(pid);
            const auto  dr  = sys.adjust_direction(pos - r_i);

            const auto dist = math::length(dr);
            E += this->potential_.potential(pid, dist);
        }
        return E;
    }

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        // update system parameters such as temperature, ionic-strength, etc.,
        this->potential_.update(sys);
        return;
    }

    void update(const system_type& sys) override
    {
        // update system parameters such as temperature, ionic-strength, etc.
        this->potential_.update(sys);
        return;
    }

    // no margin to be updated exists.
    void update_margin(const real_type, const system_type&) override {return;}

    std::string name() const override
    {return "PositionRestraint:"_s + potential_.name();}

  private:

    shape_type     shape_;
    potential_type potential_;
};

} // mjolnir
#endif//MJOLNIR_BOX_INTEARACTION_BASE
