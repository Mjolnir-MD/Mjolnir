#ifndef MJOLNIR_INPUT_READ_SIMULATOR_HPP
#define MJOLNIR_INPUT_READ_SIMULATOR_HPP
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/MolecularDynamicsSimulator.hpp>
#include <mjolnir/core/SteepestDescentSimulator.hpp>
#include <mjolnir/core/SimulatedAnnealingSimulator.hpp>
#include <mjolnir/core/SwitchingForceFieldSimulator.hpp>
#include <mjolnir/core/EnergyCalculationSimulator.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <mjolnir/util/logger.hpp>
#include <mjolnir/input/read_system.hpp>
#include <mjolnir/input/read_forcefield.hpp>
#include <mjolnir/input/read_integrator.hpp>
#include <mjolnir/input/read_observer.hpp>
#include <mjolnir/input/read_random_number_generator.hpp>
#include <mjolnir/input/read_path.hpp>

namespace mjolnir
{

template<typename traitsT, typename integratorT>
std::unique_ptr<SimulatorBase>
read_molecular_dynamics_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
        "parallelism"_s, "seed"_s, "total_step"_s, "save_step"_s, "delta_t"_s,
        "integrator"_s, "env"_s});

    const auto tstep = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep = toml::find<std::size_t>(simulator, "save_step");
    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    // later move them, so non-const
    auto sys  = read_system    <traitsT>(root, 0);
    auto obs  = read_observer  <traitsT>(root);
    auto ff   = read_forcefield<traitsT>(root, 0, simulator);
    auto rng  = read_rng       <traitsT>(simulator);
    auto intg = read_integrator<integratorT>(simulator);

    return make_unique<MolecularDynamicsSimulator<traitsT, integratorT>>(
            tstep, sstep, std::move(sys), std::move(ff), std::move(intg),
            std::move(obs), std::move(rng));
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
            read_system<traitsT>(root, 0),
            read_forcefield<traitsT>(root, 0, simulator),
            read_observer<traitsT>(root));
}

template<typename traitsT, typename integratorT>
std::unique_ptr<SimulatorBase>
read_simulated_annealing_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type   = typename traitsT::real_type;

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
            "parallelism"_s, "seed"_s, "total_step"_s, "save_step"_s, "delta_t"_s,
            "integrator"_s, "schedule"_s, "each_step"_s, "env"_s});

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

    // check integrator has temperature control
    const auto& integrator     = toml::find(simulator, "integrator");
    const auto integrator_type = toml::find<std::string>(integrator, "type");
    if(integrator_type == "VelocityVerlet")
    {
        MJOLNIR_LOG_ERROR("Simulated Annealing + NVE Newtonian");
        MJOLNIR_LOG_ERROR("NVE Newtonian doesn't have temperature control.");

        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_simulated_annealing_simulator: invalid integrator: ",
            toml::find(integrator, "type"), "here", {
            "Newtonian Integrator does not controls temperature."
            "expected value is one of the following.",
            "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                      " based on the Velocity Verlet",
            "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator"
            }));
    }

    auto sys  = read_system    <traitsT>(root, 0);
    auto ff   = read_forcefield<traitsT>(root, 0, simulator);
    auto obs  = read_observer  <traitsT>(root);
    auto rng  = read_rng       <traitsT>(simulator);
    auto intg = read_integrator<integratorT>(simulator);

    if(schedule_type == "linear")
    {
        MJOLNIR_LOG_NOTICE("temparing schedule is linear.");
        auto sch = LinearScheduler<real_type>(schedule_begin, schedule_end);

        return make_unique<
            SimulatedAnnealingSimulator<traitsT, integratorT, LinearScheduler>>(
                tstep, sstep, each_step, std::move(sch),  std::move(sys),
                std::move(ff), std::move(intg), std::move(obs), std::move(rng));
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

template<typename traitsT, typename integratorT>
std::unique_ptr<SimulatorBase>
read_switching_forcefield_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
            "parallelism"_s, "seed"_s, "total_step"_s, "save_step"_s,
            "delta_t"_s, "integrator"_s, "schedule"_s, "env"_s});

    const auto tstep = toml::find<std::size_t>(simulator, "total_step");
    const auto sstep = toml::find<std::size_t>(simulator, "save_step");
    MJOLNIR_LOG_NOTICE("total step is ", tstep);
    MJOLNIR_LOG_NOTICE("save  step is ", sstep);

    // later move them, so non-const
    auto sys  = read_system  <traitsT>(root, 0);
    auto obs  = read_observer<traitsT>(root);
    auto rng  = read_rng     <traitsT>(simulator);
    auto intg = read_integrator<integratorT>(simulator);

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
    using forcefield_type = std::unique_ptr<ForceFieldBase<traitsT>>;

    const auto forcefields = toml::find<toml::array>(root, "forcefields");
    std::vector<forcefield_type> ffs; ffs.reserve(forcefields.size());
    std::map<std::string, std::size_t> ffidx;

    for(std::size_t i=0; i<forcefields.size(); ++i)
    {
        const auto name = toml::find<std::string>(forcefields.at(i), "name");
        ffidx[name] = i;
        ffs.push_back(read_forcefield<traitsT>(root, i, simulator));

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

    return make_unique<SwitchingForceFieldSimulator<traitsT, integratorT>>(
            tstep, sstep, std::move(sys), std::move(ffs), std::move(intg),
            std::move(obs), std::move(rng), std::move(ffidx), std::move(sch));
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_energy_calculation_simulator(
        const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    using real_type       = typename traitsT::real_type;
    using coordinate_type = typename traitsT::coordinate_type;

    check_keys_available(simulator, {"type"_s, "boundary_type"_s, "precision"_s,
            "file"_s, "parallelism"_s, "env"_s});

    // ------------------------------------------------------------------------
    // construct observers manually ...
    //
    // This is the only simulator that does not need to output a trajectory,
    // therefore `XXXObserver`s are not needed except EnergyObserver.

    const auto& output = toml::find(root, "files", "output");
    const auto progress_bar_enabled = toml::find_or<bool>(output, "progress_bar", true);

    const auto output_path   = read_output_path(root);
    const auto output_prefix = toml::find<std::string>(output, "prefix");
    const auto output_name   = output_path + output_prefix;
    MJOLNIR_LOG_NOTICE("output file prefix is `", output_path, output_prefix, '`');

    // push EnergyObserver to the observer.
    ObserverContainer<traitsT> obs(progress_bar_enabled);
    obs.push_back(make_unique<EnergyObserver<traitsT>>(output_name));

    // ------------------------------------------------------------------------
    // read filename, set-up Loader

    const auto input_file = toml::find<std::string>(simulator, "file");
    std::unique_ptr<LoaderBase<traitsT>> loader;

    if(input_file.substr(input_file.size()-4, 4) == ".xyz")
    {
        loader = make_unique<XYZLoader<traitsT>>(input_file);
    }
    else if(input_file.substr(input_file.size()-4, 4) == ".dcd")
    {
        loader = make_unique<DCDLoader<traitsT>>(input_file);
    }
    else if(input_file.substr(input_file.size()-4, 4) == ".trr")
    {
        loader = make_unique<TRRLoader<traitsT>>(input_file);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_energy_calculation_simulator: invalid file extension: ",
            toml::find(simulator, "file"), "here", {
            "expected filetype is one of the following:",
            "- \"xyz\" : simple ascii file format.",
            "- \"dcd\" : binary file format.",
            "- \"trr\" : binary file format."
            }));
    }
    loader->initialize();

    // ------------------------------------------------------------------------
    // read [[systems]] manualy ...
    //
    // XXX We don't need the initial configuration, but still, the name and the
    //     group names are required to calculate energy correctly.
    //
    // This is the only simulator that loads the current state from a different
    // file. So here `positions` are not required.

    const auto& systems = toml::find(root, "systems").as_array();
    if(systems.size() != 1)
    {
        throw_exception<std::out_of_range>("[error] mjolnir::read_system: "
            "invalid number of system definitions: ", systems.size());
    }

    MJOLNIR_LOG_NOTICE("reading a system ...");
    const auto system = read_table_from_file(systems.at(0), "systems");

    const auto& boundary = toml::find(system, "boundary_shape");
    const auto& particles = toml::find(system, "particles").as_array();
    System<traitsT> sys(particles.size(), read_boundary<traitsT>(boundary));

    if(particles.size() != loader->num_particles())
    {
        throw_exception<std::runtime_error>("[error] "
            "mjolnir::read_energy_calculation_simulator: invalid number of "
            "particles: [[systems]] table has ", particles.size(),
            " elements but ", input_file, " has ", sys.size(), " particles.",
            toml::format_error("", toml::find(system, "particles"), "here"));
    }

    // clear all parameters once
    for(std::size_t i=0; i<sys.size(); ++i)
    {
        const auto& p = particles.at(i);
        // position(pos), velocity(vel) are ignored here
        check_keys_available(p, {"m"_s, "mass"_s, "name"_s, "group"_s});

        if(p.as_table().count("m") == 1)
        {
            sys.mass(i) = toml::find<real_type>(p, "m");
        }
        else
        {
            sys.mass(i) = toml::find<real_type>(p, "mass");
        }
        sys.rmass(i)    = real_type(1) / sys.mass(i);
        sys.position(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        sys.velocity(i) = math::make_coordinate<coordinate_type>(0, 0, 0);
        sys.force(i)    = math::make_coordinate<coordinate_type>(0, 0, 0);
        sys.name(i)     = toml::find_or<std::string>(p, "name",  "X"   );
        sys.group(i)    = toml::find_or<std::string>(p, "group", "NONE");
    }

    for(const auto& attr : toml::find<toml::table>(system, "attributes"))
    {
        const real_type attribute = toml::get<real_type>(attr.second);
        sys.attribute(attr.first) = attribute;
        MJOLNIR_LOG_INFO("attribute.", attr.first, " = ", attribute);
    }

    auto ff = read_forcefield<traitsT>(root, 0, simulator);

    return make_unique<EnergyCalculationSimulator<traitsT>>(loader->num_frames(),
            std::move(loader), std::move(sys), std::move(ff), std::move(obs));
}


template<typename traitsT, typename integratorT>
std::unique_ptr<SimulatorBase>
read_simulator(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();

    const auto type = toml::find<std::string>(simulator, "type");
    if(type == "MolecularDynamics")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is MolecularDynamics.");
        return read_molecular_dynamics_simulator<traitsT, integratorT>(root, simulator);
    }
    else if(type == "SteepestDescent")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is SteepestDescent.");

        // It does not do "time integration", so integratorT is not needed.
        return read_steepest_descent_simulator<traitsT>(root, simulator);
    }
    else if(type == "SimulatedAnnealing")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is SimulatedAnnealing.");
        return read_simulated_annealing_simulator<traitsT, integratorT>(root, simulator);
    }
    else if(type == "SwitchingForceField")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is SwitchingForceField.");
        return read_switching_forcefield_simulator<traitsT, integratorT>(root, simulator);
    }
    else if(type == "EnergyCalculation")
    {
        MJOLNIR_LOG_NOTICE("Simulator type is EnergyCalculation.");

        // It does not do "time integration", so integratorT is not needed.
        return read_energy_calculation_simulator<traitsT>(root, simulator);
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
            "- \"SwitchingForceField\": switch forcefield while running simulation",
            "- \"EnergyCalculation\"  : calculate energy based on a trajectory file"
            }));
    }
}

template<typename traitsT>
std::unique_ptr<SimulatorBase>
read_integrator_type(const toml::value& root, const toml::value& simulator)
{
    MJOLNIR_GET_DEFAULT_LOGGER();
    MJOLNIR_LOG_FUNCTION();
    // just determine types, not read the whole integrator.

    // There are some Simulators that does not require time integration, such as
    // SteepestDescentSimulator. In that case, any `integrator.type` is defined.
    // So if any integrator is defined, it uses VelocityVerletIntegrator for a
    // workaround.
    //     When the read_integrator called later, read_integrator function
    // throws an exception which states "no integrator defined". Simulators that
    // does not require integrator does not call read_integrator, so the error
    // can be detected.
    if(simulator.as_table().count("integrator") == 0)
    {
        using integratorT = VelocityVerletIntegrator<traitsT>;
        return read_simulator<traitsT, integratorT>(root, simulator);
    }

    // if simulator.integrator is defined, simulator.integrator.type should be
    // defined.

    const auto integ = toml::find<std::string>(simulator, "integrator", "type");
    if(integ == "VelocityVerlet")
    {
        MJOLNIR_LOG_NOTICE("Integrator is VelocityVerlet.");
        using integratorT = VelocityVerletIntegrator<traitsT>;
        return read_simulator<traitsT, integratorT>(root, simulator);
    }
    else if(integ == "UnderdampedLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is UnderdampedLangevin.");
        using integratorT = UnderdampedLangevinIntegrator<traitsT>;
        return read_simulator<traitsT, integratorT>(root, simulator);
    }
    else if(integ == "BAOABLangevin")
    {
        MJOLNIR_LOG_NOTICE("Integrator is BAOABLangevin.");
        using integratorT = BAOABLangevinIntegrator<traitsT>;
        return read_simulator<traitsT, integratorT>(root, simulator);
    }
    else
    {
        throw_exception<std::runtime_error>(toml::format_error("[error] "
            "mjolnir::read_integrator_kind: invalid integrator: ",
            toml::find(simulator, "integrator", "type"), "here", {
            "expected value is one of the following.",
            "- \"VelocityVerlet\"     : simple and standard Velocity Verlet integrator.",
            "- \"UnderdampedLangevin\": simple Underdamped Langevin Integrator"
                                      " based on the Velocity Verlet",
            "- \"BAOABLangevin\"      : well-known BAOAB Langevin Integrator"
            }));
    }
}

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
// ----------------------------------------------------------------------------
// read_molecular_dynamics_simulator

extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_molecular_dynamics_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_simulated_annealing_simulator

extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulated_annealing_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_steepest_descent_simulator

extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_steepest_descent_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_switching_forcefield_simulator

extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_switching_forcefield_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_simulator
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , VelocityVerletIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, VelocityVerletIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , UnderdampedLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, UnderdampedLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<double, UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  UnlimitedBoundary>       , BAOABLangevinIntegrator<SimulatorTraits<float,  UnlimitedBoundary>       >>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<double, CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<double, CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_simulator<SimulatorTraits<float,  CuboidalPeriodicBoundary>, BAOABLangevinIntegrator<SimulatorTraits<float,  CuboidalPeriodicBoundary>>>(const toml::value&, const toml::value&);

// ----------------------------------------------------------------------------
// read_integrator_type

extern template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<double, UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<float,  UnlimitedBoundary>       >(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<double, CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);
extern template std::unique_ptr<SimulatorBase> read_integrator_type<SimulatorTraits<float,  CuboidalPeriodicBoundary>>(const toml::value&, const toml::value&);

} // mjolnir
#endif

#endif// MJOLNIR_READ_SIMULATOR
