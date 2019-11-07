#ifndef MJOLNIR_INPUT_READ_SIMULATOR_HPP
#define MJOLNIR_INPUT_READ_SIMULATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/MolecularDynamicsSimulator.hpp>
#include <mjolnir/core/SteepestDescentSimulator.hpp>
#include <mjolnir/core/SimulatedAnnealingSimulator.hpp>
#include <mjolnir/core/SwitchingForceFieldSimulator.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_system.hpp>
#include <mjolnir/input/read_forcefield.hpp>
#include <mjolnir/input/read_integrator.hpp>
#include <mjolnir/input/read_observer.hpp>
#include <mjolnir/input/read_path.hpp>

namespace mjolnir
{

template<typename traitsT>
RandomNumberGenerator<traitsT> read_rng(const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    std::uint32_t seed = 0;
    if(simulator.as_table().count("integrator") != 0)
    {
        const auto& integrator = toml::find(simulator, "integrator");
        if(integrator.as_table().count("seed") != 0)
        {
            MJOLNIR_LOG_WARN("deprecated: put `seed` under [simulator] table.");
            MJOLNIR_LOG_WARN("deprecated: ```toml");
            MJOLNIR_LOG_WARN("deprecated: [simulator]");
            MJOLNIR_LOG_WARN("deprecated: seed = 12345");
            MJOLNIR_LOG_WARN("deprecated: ```");
            seed = toml::find<std::uint32_t>(integrator, "seed");
        }
        else
        {
            seed = toml::find<std::uint32_t>(simulator, "seed");
        }
    }
    else
    {
        seed = toml::find<std::uint32_t>(simulator, "seed");
    }
    MJOLNIR_LOG_NOTICE("seed is ", seed);
    return RandomNumberGenerator<traitsT>(seed);
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_molecular_dynamics_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
        "parallelism"_s, "seed"_s, "total_step"_s, "save_step"_s, "delta_t"_s,
        "integrator"_s});

    const auto tstep = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep = toml::find<std::size_t>(simulator, "save_step");
    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    // later move them, so non-const
    auto sys = read_system    <traitsT>(root, 0);
    auto obs = read_observer  <traitsT>(root);
    auto ff  = read_forcefield<traitsT>(root, 0);
    auto rng = read_rng       <traitsT>(simulator);

    const auto& integrator     = toml::find(simulator, "integrator");
    const auto integrator_type = toml::find<std::string>(integrator, "type");

    if(integrator_type == "VelocityVerlet")
    {
        MJOLNIR_LOG_NOTICE("Integrator is VelocityVerlet.");
        using integrator_t = VelocityVerletIntegrator<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;

        auto intg = read_velocity_verlet_integrator<traitsT>(simulator);

        return make_unique<simulator_t>(tstep, sstep, std::move(sys),
                std::move(ff), std::move(intg), std::move(obs), std::move(rng));
    }
    else if(integrator_type == "UnderdampedLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is Underdamped Langevin.");
        using integrator_t = UnderdampedLangevinIntegrator<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;

        auto intg = read_underdamped_langevin_integrator<traitsT>(simulator);

        return make_unique<simulator_t>(tstep, sstep, std::move(sys),
                std::move(ff), std::move(intg), std::move(obs), std::move(rng));
    }
    else if(integrator_type == "BAOABLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is BAOAB Langevin.");
        using integrator_t = BAOABLangevinIntegrator<traitsT>;
        using simulator_t  = MolecularDynamicsSimulator<traitsT, integrator_t>;

        auto intg = read_BAOAB_langevin_integrator<traitsT>(simulator);

        return make_unique<simulator_t>(tstep, sstep, std::move(sys),
                std::move(ff), std::move(intg), std::move(obs), std::move(rng));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_molecular_dynamics_simulator: invalid integrator: ",
            toml::find(integrator, "type"), "here", {
            "expected value is one of the following.",
            "- \"VelocityVerlet\"     : simple and standard Velocity Verlet integrator.",
            "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                      " based on the Velocity Verlet",
            "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_steepest_descent_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type      = typename traitsT::real_type;
    using simulator_type = SteepestDescentSimulator<traitsT>;

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
        "parallelism"_s, "step_limit"_s, "save_step"_s, "delta"_s, "threshold"_s});

    const auto step_lim  = toml::find<std::size_t>(simulator, "step_limit");
    const auto save_step = toml::find<std::size_t>(simulator, "save_step");
    const auto delta     = toml::find<real_type  >(simulator, "delta");
    const auto threshold = toml::find<real_type  >(simulator, "threshold");

    MJOLNIR_LOG_NOTICE("step_limit is ", step_lim);
    MJOLNIR_LOG_NOTICE("save_step  is ", save_step);
    MJOLNIR_LOG_NOTICE("delta      is ", delta);
    MJOLNIR_LOG_NOTICE("threshold  is ", threshold);

    return make_unique<simulator_type>(delta, threshold, step_lim, save_step,
            read_system<traitsT>(root, 0), read_forcefield<traitsT>(root, 0),
            read_observer<traitsT>(root));
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulated_annealing_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type   = typename traitsT::real_type;

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
            "parallelism"_s, "seed"_s, "total_step"_s, "save_step"_s,
            "delta_t"_s, "integrator"_s, "schedule"_s, "each_step"_s});

    const auto tstep = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep = toml::find<std::size_t>(simulator, "save_step");

    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    const auto schedule       = toml::find(simulator, "schedule");
    const auto schedule_type  = toml::find<std::string>(schedule, "type");
    const auto schedule_begin = toml::find<real_type>  (schedule, "begin");
    const auto schedule_end   = toml::find<real_type>  (schedule, "end");
    const auto each_step      = toml::find<std::size_t>(simulator, "each_step");

    check_keys_available(schedule, {"type"_s, "begin"_s, "end"_s, "each_step"_s});

    MJOLNIR_LOG_NOTICE("temperature from ", schedule_begin);
    MJOLNIR_LOG_NOTICE("temperature to   ", schedule_end);
    MJOLNIR_LOG_INFO("update temperature for each ", each_step, " steps");

    const auto& integrator     = toml::find(simulator, "integrator");
    const auto integrator_type = toml::find<std::string>(integrator, "type");

    auto sys = read_system    <traitsT>(root, 0);
    auto ff  = read_forcefield<traitsT>(root, 0);
    auto obs = read_observer  <traitsT>(root);
    auto rng = read_rng       <traitsT>(simulator);

    if(schedule_type == "linear")
    {
        MJOLNIR_LOG_NOTICE("temparing schedule is linear.");
        auto sch = LinearScheduler<real_type>(schedule_begin, schedule_end);

        if(integrator_type == "VelocityVerlet")
        {
            MJOLNIR_LOG_ERROR("Simulated Annealing + NVE Newtonian");
            MJOLNIR_LOG_ERROR("NVE Newtonian doesn't have temperature control.");

            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_simulated_annealing_simulator: invalid integrator: ",
                toml::find(integrator, "type"), "here", {
                "Newtonian Integrator does not controls temperature. "
                "expected value is one of the following.",
                "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                          " based on the Velocity Verlet",
                "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator"
                }));
        }
        else if(integrator_type == "UnderdampedLangevin")
        {
            MJOLNIR_LOG_NOTICE("Integrator is Underdamped Langevin.");
            using integrator_t = UnderdampedLangevinIntegrator<traitsT>;
            using simulator_t  = SimulatedAnnealingSimulator<
                traitsT, integrator_t, LinearScheduler>;

            auto intg = read_underdamped_langevin_integrator<traitsT>(simulator);

            return make_unique<simulator_t>(tstep, sstep, each_step,
                    std::move(sch),  std::move(sys), std::move(ff),
                    std::move(intg), std::move(obs), std::move(rng));
        }
        else if(integrator_type == "BAOABLangevin")
        {
            MJOLNIR_LOG_NOTICE("Integrator is BAOAB Langevin.");
            using integrator_t = BAOABLangevinIntegrator<traitsT>;
            using simulator_t  = SimulatedAnnealingSimulator<
                traitsT, integrator_t, LinearScheduler>;

            auto intg = read_BAOAB_langevin_integrator<traitsT>(simulator);

            return make_unique<simulator_t>(tstep, sstep, each_step,
                    std::move(sch),  std::move(sys), std::move(ff),
                    std::move(intg), std::move(obs), std::move(rng));
        }
        else
        {
            throw_exception<std::runtime_error>(toml::format_error("[error] "
                "mjolnir::read_simulated_annealing_simulator: invalid integrator: ",
                toml::find(simulator, "integrator"), "here", {
                "expected value is one of the following.",
                "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                          " based on the Velocity Verlet"
                }));
        }
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_simulated_annealing_simulator: invalid schedule.type",
            toml::find(schedule, "type"), "here", {
            "expected value is one of the following.",
            "- \"linear\"     : simple linear temperature scheduling"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_switching_forcefield_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
            "parallelism"_s, "seed"_s, "total_step"_s, "save_step"_s,
            "delta_t"_s, "integrator"_s, "schedule"_s});

    const auto tstep = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep = toml::find<std::size_t>(simulator, "save_step");
    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    // later move them, so non-const
    auto sys = read_system  <traitsT>(root, 0);
    auto obs = read_observer<traitsT>(root);
    auto rng = read_rng     <traitsT>(simulator);

    // ------------------------------------------------------------------------
    // read schedule
    const auto& schedule = toml::find(simulator, "schedule").as_array();
    std::vector<std::pair<std::size_t, std::string>> sch(schedule.size());
    std::transform(schedule.begin(), schedule.end(), sch.begin(),
            [](const toml::value& v){
                return std::make_pair(toml::find<std::size_t>(v, "until"),
                                      toml::find<std::string>(v, "forcefield"));
            });

    // ------------------------------------------------------------------------
    // read all forcefields and its names
    using forcefield_type = ForceField<traitsT>;

    const auto forcefields = toml::find<toml::array>(root, "forcefields");
    std::vector<forcefield_type> ffs; ffs.reserve(forcefields.size());
    std::map<std::string, std::size_t> ffidx;

    for(std::size_t i=0; i<forcefields.size(); ++i)
    {
        const auto name = toml::find<std::string>(forcefields.at(i), "name");
        ffidx[name] = i;
        ffs.push_back(read_forcefield<traitsT>(root, i));

        MJOLNIR_LOG_NOTICE(i, "-th forcefield \"", name, "\" found.");
    }
    MJOLNIR_LOG_NOTICE("forcefields = ", ffidx);

    // ------------------------------------------------------------------------
    // check consistency
    for(const auto& v : schedule)
    {
        const auto& ff_name = toml::find<std::string>(v, "forcefield");
        if(ffidx.count(ff_name) == 0)
        {
            MJOLNIR_LOG_ERROR("forcefield \"", ff_name, "\" does not exist.");
            throw std::runtime_error(toml::format_error("[error] "
                    "forcefield is not defined.", toml::find(v, "forcefield"),
                    "this is not defined in [[forcefields]]"));
        }
    }

    // ------------------------------------------------------------------------
    // read integrator and then return simulator

    const auto& integrator     = toml::find(simulator, "integrator");
    const auto integrator_type = toml::find<std::string>(integrator, "type");

    if(integrator_type == "VelocityVerlet")
    {
        MJOLNIR_LOG_NOTICE("Integrator is VelocityVerlet.");
        using integrator_t = VelocityVerletIntegrator<traitsT>;
        using simulator_t  = SwitchingForceFieldSimulator<traitsT, integrator_t>;

        auto intg = read_velocity_verlet_integrator<traitsT>(simulator);

        return make_unique<simulator_t>(tstep, sstep, std::move(sys),
                std::move(ffs), std::move(intg), std::move(obs),
                std::move(rng), std::move(ffidx), std::move(sch));
    }
    else if(integrator_type == "UnderdampedLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is Underdamped Langevin.");
        using integrator_t = UnderdampedLangevinIntegrator<traitsT>;
        using simulator_t  = SwitchingForceFieldSimulator<traitsT, integrator_t>;

        auto intg = read_underdamped_langevin_integrator<traitsT>(simulator);

        return make_unique<simulator_t>(tstep, sstep, std::move(sys),
                std::move(ffs), std::move(intg), std::move(obs),
                std::move(rng), std::move(ffidx), std::move(sch));
    }
    else if(integrator_type == "BAOABLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is BAOAB Langevin.");
        using integrator_t = BAOABLangevinIntegrator<traitsT>;
        using simulator_t  = SwitchingForceFieldSimulator<traitsT, integrator_t>;

        auto intg = read_BAOAB_langevin_integrator<traitsT>(simulator);

        return make_unique<simulator_t>(tstep, sstep, std::move(sys),
                std::move(ffs), std::move(intg), std::move(obs),
                std::move(rng), std::move(ffidx), std::move(sch));
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_switching_forcefield_simulator: invalid integrator: ",
            toml::find(integrator, "type"), "here", {
            "expected value is one of the following.",
            "- \"VelocityVerlet\"     : simple and standard Velocity Verlet integrator.",
            "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                      " based on the Velocity Verlet",
            "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_simulator(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto type = toml::find<std::string>(simulator, "type");
    if(type == "MolecularDynamics")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is MolecularDynamics.");
        return read_molecular_dynamics_simulator<traitsT>(root, simulator);
    }
    else if(type == "SteepestDescent")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is SteepestDescent.");
        return read_steepest_descent_simulator<traitsT>(root, simulator);
    }
    else if(type == "SimulatedAnnealing")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is SimulatedAnnealing.");
        return read_simulated_annealing_simulator<traitsT>(root, simulator);
    }
    else if(type == "SwitchingForceField")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is SwitchingForceField.");
        return read_switching_forcefield_simulator<traitsT>(root, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_simulator: invalid type",
            toml::find(simulator, "type"), "here", {
            "expected value is one of the following.",
            "- \"MolecularDynamcis\"  : standard MD simulation",
            "- \"SteepestDescent\"    : energy minimization by gradient method",
            "- \"SimulatedAnnealing\" : energy minimization by Annealing",
            "- \"SwitchingForceField\": switch forcefield while running simulation"
            }));
    }
}

} // mjolnir


#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

extern template RandomNumberGenerator<SimulatorTraits<double, UnlimitedBoundary>       > read_rng(const toml::value&);
extern template RandomNumberGenerator<SimulatorTraits<float,  UnlimitedBoundary>       > read_rng(const toml::value&);
extern template RandomNumberGenerator<SimulatorTraits<double, CuboidalPeriodicBoundary>> read_rng(const toml::value&);
extern template RandomNumberGenerator<SimulatorTraits<float,  CuboidalPeriodicBoundary>> read_rng(const toml::value&);
} // mjolnir
#endif

#endif// MJOLNIR_READ_SIMULATOR
