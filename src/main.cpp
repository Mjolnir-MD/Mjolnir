#include <mjolnir/core/Simulator.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/io/toml.hpp>
typedef mjolnir::DefaultTraits traits;

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input.toml>" << std::endl;
        return 1;
    }

    std::ifstream ifs(argv[1]);
    if(!ifs.good())
    {
        std::cerr << "file open error: " << argv[1] << std::endl;
        return 0;
    }

    // read input file
    auto input = toml::parse(ifs);

    // read [simulator] block
    auto sim    = toml::get<toml::Table>(input.at("simulator"));
    auto pcon   = mjolnir::read_particles<traits>(sim);
    auto rng    = mjolnir::read_random_number_generator<traits>(sim);
    auto integr = mjolnir::read_time_integrator<traits>(sim, rng);
    const std::size_t total_step = toml::get<toml::Integer>(sim.at("total_step"));
    const std::size_t save_step  = toml::get<toml::Integer>(sim.at("save_step"));

    const std::string filename = toml::get<toml::String>(sim.at("file_name"));
    const std::string trajname = filename + ".xyz";
    const std::string ene_name = filename + ".ene";
    mjolnir::Observer<traits> obs(trajname, ene_name);

//     auto forcefield = toml::get<toml::Table>(input.at("forcefield"));
    mjolnir::ForceField<traits> ff = mjolnir::read_force_field<traits>(input);

    mjolnir::Simulator<traits> simulator(
            std::move(pcon), std::move(ff), std::move(integr));

    // run md
    std::cerr << "start running simulation" << std::endl;
    simulator.initialize();
    for(std::size_t i=0; i<total_step; ++i)
    {
        if(i % save_step == 0)
        {
            obs.output_coordinate(simulator);
            obs.output_energy(simulator);
        }
        simulator.step();
    }

    // output last state
    obs.output_coordinate(simulator);
    obs.output_energy(simulator);

    return 0;
}
