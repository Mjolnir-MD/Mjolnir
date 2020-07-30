#ifndef MJOLNIR_INTERACTION_EXTERNAL_RECTANGULAR_BOX_INTERACTION_HPP
#define MJOLNIR_INTERACTION_EXTERNAL_RECTANGULAR_BOX_INTERACTION_HPP
#include <mjolnir/core/ExternalForceInteractionBase.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/math/math.hpp>
#include <mjolnir/util/format_nth.hpp>

namespace mjolnir
{

/*! @brief Interaction between particle and a box. */
template<typename traitsT, typename potentialT>
class RectangularBoxInteraction final
    : public ExternalForceInteractionBase<traitsT>
{
  public:
    using traits_type     = traitsT;
    using potential_type  = potentialT;
    using base_type       = ExternalForceInteractionBase<traitsT>;
    using real_type       = typename base_type::real_type;
    using coordinate_type = typename base_type::coordinate_type;
    using system_type     = typename base_type::system_type;
    using boundary_type   = typename base_type::boundary_type;

  public:

    RectangularBoxInteraction(
        const coordinate_type& lower, const coordinate_type& upper,
        const real_type       margin, potential_type&& pot)
        : cutoff_(-1), margin_(margin), current_margin_(-1),
          lower_(lower), upper_(upper), neighbors_{}, potential_(std::move(pot))
    {
        MJOLNIR_GET_DEFAULT_LOGGER();

        if(!is_unlimited_boundary<boundary_type>::value)
        {
            const auto msg = "RectangularBox cannot be used under the periodic "
                             "boundary condition. Use unlimited boundary.";
            MJOLNIR_LOG_ERROR(msg);
            throw std::runtime_error(msg);
        }
        assert(math::X(lower_) < math::X(upper_));
        assert(math::Y(lower_) < math::Y(upper_));
        assert(math::Z(lower_) < math::Z(upper_));
    }
    ~RectangularBoxInteraction() override {}

    // calculate force, update spatial partition (reduce margin) inside.
    void      calc_force (system_type&)           const noexcept override;
    real_type calc_energy(system_type const&)     const noexcept override;
    real_type calc_force_and_energy(system_type&) const noexcept override;

    /*! @brief initialize spatial partition (e.g. CellList)                   *
     *  @details before calling `calc_(force|energy)`, this should be called. */
    void initialize(const system_type& sys) override
    {
        this->potential_.update(sys); // update system parameters
        this->cutoff_ = this->potential_.max_cutoff_length();
        this->make(sys);              // update neighbor list
        return;
    }

    /*! @brief update parameters (e.g. temperature, ionic strength, ...)  *
     *  @details A method that change system parameters (e.g. Annealing), *
     *           the method is bound to call this function after changing *
     *           parameters.                                              */
    void update(const system_type& sys) override
    {
        this->potential_.update(sys); // update system parameters
        this->make(sys);              // update neighbor list
    }

    void reduce_margin(const real_type dmargin, const system_type& sys) override
    {
        this->current_margin_ -= dmargin;
        if(this->current_margin_ < 0)
        {
            this->make(sys);
        }
        return;
    }
    void scale_margin(const real_type scale, const system_type& sys) override
    {
        this->current_margin_ = (cutoff_ + current_margin_) * scale - cutoff_;
        if(this->current_margin_ < 0)
        {
            this->make(sys);
        }
        return;
    }

    std::string name() const override
    {return "RectangularBox:"_s + potential_.name();}

    base_type* clone() const override
    {
        return new RectangularBoxInteraction(lower_, upper_, margin_,
                                             potential_type(potential_));
    }

    // for tests
    coordinate_type const& upper() const noexcept {return upper_;}
    coordinate_type const& lower() const noexcept {return lower_;}
    potential_type  const& potential() const noexcept {return potential_;}

  private:

    // construct a neighbor-list.
    // Also, check if all the particles are inside of the box.
    void make(const system_type& sys)
    {
        MJOLNIR_GET_DEFAULT_LOGGER();
        const real_type threshold = this->cutoff_ * (1 + this->margin_);

        this->neighbors_.clear();
        for(std::size_t i : this->potential_.participants())
        {
            const auto& pos  = sys.position(i);
            // Assuming that PBC and Box interaction will not be used together
            const auto  dr_u = this->upper_ - pos;
            const auto  dr_l = pos - this->lower_;

            // check if all the particles are inside of the box.
            if(math::X(dr_u) < 0 || math::Y(dr_u) < 0 || math::Z(dr_u) < 0)
            {
                MJOLNIR_LOG_ERROR(format_nth(i), " particle, ", pos,
                                  " exceeds the upper bound, ", this->upper_);
                throw_exception<std::runtime_error>(format_nth(i), " particle ",
                        pos, " exceeds the upper bound, ", this->upper_);
            }
            if(math::X(dr_l) < 0 || math::Y(dr_l) < 0 || math::Z(dr_l) < 0)
            {
                MJOLNIR_LOG_ERROR(format_nth(i), " particle, ", pos,
                                  " exceeds the lower bound, ", this->lower_);
                throw_exception<std::runtime_error>(format_nth(i), " particle ",
                        pos, " exceeds the lower bound, ", this->lower_);

            }

            // assign particles that are within the cutoff range.
            if(math::X(dr_u) <= threshold || math::X(dr_l) <= threshold ||
               math::Y(dr_u) <= threshold || math::Y(dr_l) <= threshold ||
               math::Z(dr_u) <= threshold || math::Z(dr_l) <= threshold)
            {
                neighbors_.push_back(i);
            }
        }
        this->current_margin_ = this->cutoff_ * this->margin_;
        return;
    }

  private:

    real_type                cutoff_, margin_, current_margin_;
    coordinate_type          lower_, upper_;
    std::vector<std::size_t> neighbors_;
    potential_type           potential_;

// #ifdef MJOLNIR_WITH_OPENMP
//     // OpenMP implementation uses its own implementation to run it in parallel.
//     // So this implementation should not be instanciated with OpenMP Traits.
//     static_assert(!is_openmp_simulator_traits<traits_type>::value,
//                   "this is the default implementation, not for OpenMP");
// #endif
};

template<typename traitsT, typename potT>
void RectangularBoxInteraction<traitsT, potT>::calc_force(
        system_type& sys) const noexcept
{
    // Here we assume that all the particles are inside of the box.
    // Also we assume that no boundary condition is applied.
    for(const std::size_t i : this->neighbors_)
    {
        const auto& pos = sys.position(i);
        const auto dr_u = this->upper_ - pos;
        const auto dr_l = pos - this->lower_;

        const auto dx_u = this->potential_.derivative(i, math::X(dr_u));
        const auto dx_l = this->potential_.derivative(i, math::X(dr_l));
        const auto dy_u = this->potential_.derivative(i, math::Y(dr_u));
        const auto dy_l = this->potential_.derivative(i, math::Y(dr_l));
        const auto dz_u = this->potential_.derivative(i, math::Z(dr_u));
        const auto dz_l = this->potential_.derivative(i, math::Z(dr_l));

        sys.force(i) += math::make_coordinate<coordinate_type>(
                dx_u - dx_l, dy_u - dy_l, dz_u - dz_l);
    }
    return ;
}

template<typename traitsT, typename potT>
typename RectangularBoxInteraction<traitsT, potT>::real_type
RectangularBoxInteraction<traitsT, potT>::calc_energy(
        const system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const std::size_t i : this->neighbors_)
    {
        const auto& pos = sys.position(i);
        const auto dr_u = this->upper_ - pos;
        const auto dr_l = pos - this->lower_;

        E += this->potential_.potential(i, math::X(dr_u)) +
             this->potential_.potential(i, math::X(dr_l)) +
             this->potential_.potential(i, math::Y(dr_u)) +
             this->potential_.potential(i, math::Y(dr_l)) +
             this->potential_.potential(i, math::Z(dr_u)) +
             this->potential_.potential(i, math::Z(dr_l));
    }
    return E;
}

template<typename traitsT, typename potT>
typename RectangularBoxInteraction<traitsT, potT>::real_type
RectangularBoxInteraction<traitsT, potT>::calc_force_and_energy(
        system_type& sys) const noexcept
{
    real_type E = 0.0;
    for(const std::size_t i : this->neighbors_)
    {
        const auto& pos = sys.position(i);
        const auto dr_u = this->upper_ - pos;
        const auto dr_l = pos - this->lower_;

        const auto dx_u = this->potential_.derivative(i, math::X(dr_u));
        const auto dx_l = this->potential_.derivative(i, math::X(dr_l));
        const auto dy_u = this->potential_.derivative(i, math::Y(dr_u));
        const auto dy_l = this->potential_.derivative(i, math::Y(dr_l));
        const auto dz_u = this->potential_.derivative(i, math::Z(dr_u));
        const auto dz_l = this->potential_.derivative(i, math::Z(dr_l));

        sys.force(i) += math::make_coordinate<coordinate_type>(
                dx_u - dx_l, dy_u - dy_l, dz_u - dz_l);

        E += this->potential_.potential(i, math::X(dr_u)) +
             this->potential_.potential(i, math::X(dr_l)) +
             this->potential_.potential(i, math::Y(dr_u)) +
             this->potential_.potential(i, math::Y(dr_l)) +
             this->potential_.potential(i, math::Z(dr_u)) +
             this->potential_.potential(i, math::Z(dr_l));
    }
    return E;
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/forcefield/external/ExcludedVolumeWallPotential.hpp>
#include <mjolnir/forcefield/external/LennardJonesWallPotential.hpp>

namespace mjolnir
{

extern template class RectangularBoxInteraction<SimulatorTraits<double, UnlimitedBoundary>,        ExcludedVolumeWallPotential<double>>;
extern template class RectangularBoxInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        ExcludedVolumeWallPotential<float >>;
extern template class RectangularBoxInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<double>>;
extern template class RectangularBoxInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, ExcludedVolumeWallPotential<float >>;

extern template class RectangularBoxInteraction<SimulatorTraits<double, UnlimitedBoundary>,        LennardJonesWallPotential<double>>;
extern template class RectangularBoxInteraction<SimulatorTraits<float,  UnlimitedBoundary>,        LennardJonesWallPotential<float >>;
extern template class RectangularBoxInteraction<SimulatorTraits<double, CuboidalPeriodicBoundary>, LennardJonesWallPotential<double>>;
extern template class RectangularBoxInteraction<SimulatorTraits<float,  CuboidalPeriodicBoundary>, LennardJonesWallPotential<float >>;

} // mjolnir
#endif // MJOLNIR_SEPARATE_BUILD

#endif//MJOLNIR_BOX_INTEARACTION_BASE
