#include <mjolnir/input/read_input_file.hpp>

namespace mjolnir
{

// This potential function is the same to the case of BAOAB integrator test.
// Fundamental integration algorithm is same to BAOAB integrator, so the result 
// should be same to BAOAB's one.
template<typename realT>
class TestPotential
{
  public:
    using real_type = realT;

  public:

    explicit TestPotential(std::vector<std::size_t>&& ps)
        : participants_(std::move(ps))
    {}
    ~TestPotential() = default;

    real_type potential(const std::size_t, const real_type r) const noexcept
    {
        return 0.25 * (r * r * r * r) + std::sin(1 + 5 * r);
    }

    real_type derivative(const std::size_t, const real_type r) const noexcept
    {
        return (r * r * r) + 5 * std::cos(1 + 5 * r);
    }

    real_type max_cutoff_length() const noexcept
    {
        return std::numeric_limits<real_type>::max();
    }

    std::vector<std::size_t> participants() const
    {
        return this->participants_;
    }

    // nothing to do when system parameters change.
    template<typename T>
    void update(const System<T>&) const noexcept {return;}

    const char* name() const noexcept {return "Test";}

  private:

    std::vector<std::size_t> participants_;
};
} // mjolnir

// this function tries to inject a TestPotential and return true if succeed.
template<typename traitsT>
bool inject_test_potential(std::unique_ptr<mjolnir::SimulatorBase>& sim_base)
{
    using traits_type     = traitsT;
    using integrator_type = mjolnir::gBAOABLangevinIntegrator<traits_type>;
    using simulator_type  = mjolnir::MolecularDynamicsSimulator<traits_type, integrator_type>;

    using shape_type       = mjolnir::AxisAlignedPlane<traits_type, mjolnir::PositiveXDirection<traits_type>>;
    using potential_type   = mjolnir::TestPotential<typename traits_type::real_type>;
    using interaction_type = mjolnir::ExternalDistanceInteraction<traits_type, potential_type, shape_type>;

    std::vector<std::size_t> ps(100);
    std::iota(ps.begin(), ps.end(), 0ul);

    simulator_type* sim = dynamic_cast<simulator_type*>(sim_base.get());
    if(!sim)
    {
        return false;
    }

    sim->forcefields().external().emplace(mjolnir::make_unique<interaction_type>(
            shape_type(0.0, 1.0), potential_type(std::move(ps))));

    return true;
}

int main()
{
    // -----------------------------------------------------------------------
    // read base setup

    auto sim_base = mjolnir::read_input_file("test_g-BAOABIntegrator.toml");

    // -----------------------------------------------------------------------
    // setup external forcefield for testing and inject it into the simulator

    // check default (sequencial) implementation
    bool injected = inject_test_potential<mjolnir::SimulatorTraits<
        double, mjolnir::CuboidalPeriodicBoundary>>(sim_base);

#ifdef MJOLNIR_WITH_OPENMP
    if(!injected)
    {
        // check OpenMP implementation
        injected = inject_test_potential<mjolnir::OpenMPSimulatorTraits<
            double, mjolnir::CuboidalPeriodicBoundary>>(sim_base);
    }
#endif

    if(!injected)
    {
        std::cerr << "-- [error] the injection of TestPotential failed\n"
                     "-- Please check the simulation setup." << std::endl;
        return 1;
    }

    // -----------------------------------------------------------------------
    // run simulation

    sim_base->initialize();
    sim_base->run();
    sim_base->finalize();




    return 0;
}
