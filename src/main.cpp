#include <mjolnir/core/Simulator.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/core/VelocityVerlet.hpp>
#include <mjolnir/core/UnderdampedLangevin.hpp>
#include <mjolnir/core/BondLengthInteraction.hpp>
#include <mjolnir/core/BondAngleInteraction.hpp>
#include <mjolnir/core/DihedralAngleInteraction.hpp>
#include <mjolnir/core/HarmonicPotential.hpp>
#include <mjolnir/core/ClementiDihedralPotential.hpp>
#include <mjolnir/core/Go1012ContactPotential.hpp>
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/ExcludedVolumePotential.hpp>
#include <mjolnir/core/LennardJonesPotential.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

#include <toml/toml.hpp>

template<typename traitsT>
mjolnir::ParticleContainer<traitsT>
read_particles(const toml::Table& sim)
{
    const std::size_t num_particles =
        toml::get<toml::Integer>(sim.at("number_of_particle"));
    const std::vector<toml::Table> particles =
        toml::get<toml::Array<toml::Table>>(sim.at("particles"));
    if(num_particles != particles.size())
        throw std::runtime_error(
                "number_of_particles is not equal to particles.size()");

    mjolnir::ParticleContainer<traitsT> pcon(particles.size());

    for(auto iter = mjolnir::make_zip(particles.cbegin(), pcon.begin());
            iter != mjolnir::make_zip(particles.cend(), pcon.end()); ++iter)
    {
        const typename traitsT::real_type mass =
            toml::get<toml::Float>(mjolnir::get<0>(iter)->at("m"));
        std::vector<typename traitsT::real_type> pos = 
            toml::get<toml::Array<toml::Float>>(mjolnir::get<0>(iter)->at("r"));
        std::vector<typename traitsT::real_type> vel =
            toml::get<toml::Array<toml::Float>>(mjolnir::get<0>(iter)->at("v"));
        std::vector<typename traitsT::real_type> acc =
            toml::get<toml::Array<toml::Float>>(mjolnir::get<0>(iter)->at("a"));

        const typename traitsT::coordinate_type position(pos.at(0), pos.at(1), pos.at(2));
        const typename traitsT::coordinate_type velocity(vel.at(0), vel.at(1), vel.at(2));
        const typename traitsT::coordinate_type acceleration(acc.at(0), acc.at(1), acc.at(2));

        *mjolnir::get<1>(iter) =
            mjolnir::make_particle(mass, position, velocity, acceleration);
    }
    return pcon;
}

template<typename traitsT>
std::shared_ptr<mjolnir::RandomNumberGenerator<traitsT>>
read_random_number_generator(const toml::Table& sim)
{
    const std::uint32_t seed = toml::get<toml::Integer>(sim.at("seed"));
    return std::make_shared<mjolnir::RandomNumberGenerator<traitsT>>(seed);
}

template<typename traitsT>
std::unique_ptr<mjolnir::Integrator<traitsT>>
read_time_integrator(const toml::Table& sim,
        const std::shared_ptr<mjolnir::RandomNumberGenerator<traitsT>>& rng)
{
    const std::string integral = toml::get<toml::String>(sim.at("time_integration"));
    if(integral == "Underdamped Langevin")
    {
        const typename traitsT::real_type delta_t =
            toml::get<toml::Float>(sim.at("delta_t"));

        const typename traitsT::real_type temperature = 
            toml::get<toml::Float>(sim.at("temperature"));

        const typename traitsT::real_type kB = 
            toml::get<toml::Float>(sim.at("kB"));

        const std::size_t num_particles =
            toml::get<toml::Integer>(sim.at("number_of_particle"));

        std::vector<typename traitsT::real_type> fric =
            toml::get<toml::Array<toml::Float>>(sim.at("friction_constant"));
        if(num_particles != fric.size())
            throw std::runtime_error("friction constant not enough");

        return mjolnir::make_unique<mjolnir::UnderdampedLangevin<traitsT>>(
            delta_t, num_particles, temperature, kB, std::move(fric), rng);
    }
    else if(integral == "NVE Newtonian")
    {
        const typename traitsT::real_type delta_t =
            toml::get<toml::Float>(sim.at("delta_t"));

        const std::size_t num_particles =
            toml::get<toml::Integer>(sim.at("number_of_particle"));

        return mjolnir::make_unique<mjolnir::VelocityVerlet<traitsT>>(
                delta_t, num_particles);
    }
    else
        throw std::invalid_argument("unknown time integral method: " + integral);
}

template<typename traitsT>
std::unique_ptr<mjolnir::LocalPotentialBase<traitsT>>
read_local_potential(const std::string& name, const toml::Table& param)
{
    if(name == "Harmonic")
    {
        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));
        const typename traitsT::real_type r0 =
            toml::get<toml::Float>(param.at("r0"));
        return mjolnir::make_unique<mjolnir::HarmonicPotential<traitsT>>(k, r0);
    }
    else if(name == "ClementiDihedral")
    {
        const typename traitsT::real_type k1 =
            toml::get<toml::Float>(param.at("k1"));
        const typename traitsT::real_type k3 =
            toml::get<toml::Float>(param.at("k3"));
        const typename traitsT::real_type phi0 =
            toml::get<toml::Float>(param.at("phi0"));
        return mjolnir::make_unique<mjolnir::ClementiDihedralPotential<traitsT>>(k1, k3, phi0);
    }
    else if(name == "Go1012Contact")
    {
        const typename traitsT::real_type k =
            toml::get<toml::Float>(param.at("k"));
        const typename traitsT::real_type r0 =
            toml::get<toml::Float>(param.at("r0"));
        return mjolnir::make_unique<mjolnir::Go1012ContactPotential<traitsT>>(k, r0);
    }
    else
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
mjolnir::LocalForceField<traitsT>
read_local_force_field(const toml::Array<toml::Table>& lffs)
{
    mjolnir::LocalForceField<traitsT> lff;
    for(auto iter = lffs.cbegin(); iter != lffs.cend(); ++iter)
    {
        const std::string potential =
            toml::get<toml::String>(iter->at("potential"));
        const toml::Array<toml::Table> params = 
            toml::get<toml::Array<toml::Table>>(iter->at("parameters"));

        const std::string interaction =
            toml::get<toml::String>(iter->at("interaction"));
        if(interaction == "BondLength")
        {
            for(auto p = params.cbegin(); p != params.cend(); ++p)
            {
                std::unique_ptr<mjolnir::LocalPotentialBase<traitsT>> pot =
                    read_local_potential<traitsT>(potential, *p);
                const toml::Array<toml::Integer> indices =
                    toml::get<toml::Array<toml::Integer>>(p->at("indices"));
                if(indices.size() != 2)
                    throw std::runtime_error("bond index must be specified");
                lff.emplace_bond(indices.at(0), indices.at(1), std::move(pot));
            }
        }
        else if(interaction == "BondAngle")
        {
            for(auto p = params.cbegin(); p != params.cend(); ++p)
            {
                std::unique_ptr<mjolnir::LocalPotentialBase<traitsT>> pot =
                    read_local_potential<traitsT>(potential, *p);
                const toml::Array<toml::Integer> indices =
                    toml::get<toml::Array<toml::Integer>>(p->at("indices"));
                if(indices.size() != 3)
                    throw std::runtime_error("bond index must be specified");
                lff.emplace_angle(indices.at(0), indices.at(1), indices.at(2),
                                  std::move(pot));
            }
        }
        else if(interaction == "DihedralAngle")
        {
            for(auto p = params.cbegin(); p != params.cend(); ++p)
            {
                std::unique_ptr<mjolnir::LocalPotentialBase<traitsT>> pot =
                    read_local_potential<traitsT>(potential, *p);
                const toml::Array<toml::Integer> indices =
                    toml::get<toml::Array<toml::Integer>>(p->at("indices"));
                if(indices.size() != 4)
                    throw std::runtime_error("bond index must be specified");
                lff.emplace_dihedral(indices.at(0), indices.at(1), indices.at(2),
                        indices.at(3), std::move(pot));
            }
        }
        else
            throw std::runtime_error("unknown interaction: " + interaction);
    }
    return lff;
}

template<typename traitsT>
std::unique_ptr<mjolnir::GlobalPotentialBase<traitsT>>
read_global_potential(const std::string& name, const toml::Table& potent)
{
    if(name == "ExcludedVolume")
    {
        toml::Array<toml::Table> param_table =
            toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
        std::vector<typename mjolnir::ExcludedVolumePotential<traitsT>::parameter_type>
            params(param_table.size());
        for(auto iter = mjolnir::make_zip(param_table.cbegin(), params.begin());
                iter != mjolnir::make_zip(param_table.cend(), params.end()); ++iter)
        {
            const typename traitsT::real_type sigma = 
                toml::get<toml::Float>(mjolnir::get<0>(iter)->at("sigma"));

            *mjolnir::get<1>(iter) = sigma;
        }
        const typename traitsT::real_type epsilon = 
            toml::get<toml::Float>(potent.at("epsilon"));

        return mjolnir::make_unique<mjolnir::ExcludedVolumePotential<traitsT>>(
                epsilon, std::move(params));
    }
    else if(name == "LennardJones")
    {
        toml::Array<toml::Table> param_table =
            toml::get<toml::Array<toml::Table>>(potent.at("parameters"));
        std::vector<typename mjolnir::LennardJonesPotential<traitsT>::parameter_type>
            params(param_table.size());

        for(auto iter = mjolnir::make_zip(param_table.cbegin(), params.begin());
                iter != mjolnir::make_zip(param_table.cend(), params.end()); ++iter)
        {
            const typename traitsT::real_type sigma = 
                toml::get<toml::Float>(mjolnir::get<0>(iter)->at("sigma"));
            const typename traitsT::real_type epsilon = 
                toml::get<toml::Float>(mjolnir::get<0>(iter)->at("epsilon"));

            *mjolnir::get<1>(iter) = std::make_pair(sigma, epsilon);
        }

        return mjolnir::make_unique<mjolnir::LennardJonesPotential<traitsT>>(
                std::move(params));
    }
    else 
        throw std::runtime_error("unknown potential: " + name);
}

template<typename traitsT>
std::unique_ptr<mjolnir::GlobalInteractionBase<traitsT>>
read_global_interaction(const std::string& name,
        const std::unique_ptr<mjolnir::GlobalPotentialBase<traitsT>>& pot,
        const toml::Table& potent)
{
    if(name == "Global")
    {
        const typename traitsT::real_type cutoff = pot->max_cutoff_length();

        auto space = mjolnir::make_unique<mjolnir::VerletList<traitsT>>(
                cutoff, cutoff * 0.5);

        std::vector<std::vector<std::size_t>> excepts;
        auto elist = toml::get<toml::Array<toml::Array<toml::Integer>>>(potent.at("excepts"));
        for(auto iter = elist.cbegin(); iter != elist.cend(); ++iter)
        {
            std::vector<std::size_t> l; l.reserve(iter->size());
            for(auto j = iter->cbegin(); j != iter->cend(); ++j)
                l.push_back(static_cast<std::size_t>(*j));
            excepts.emplace_back(l);
        }
        space->set_except(excepts);

        return mjolnir::make_unique<mjolnir::GlobalDistanceInteraction<traitsT>>(
                std::move(space));
    }
    else
        throw std::runtime_error("unknown interaction: " + name);
}

template<typename traitsT>
mjolnir::GlobalForceField<traitsT>
read_global_force_field(const toml::Array<toml::Table>& gffs)
{
    mjolnir::GlobalForceField<traitsT> gff;

    for(auto iter = gffs.cbegin(); iter != gffs.cend(); ++iter)
    {
        const std::string potential = 
            toml::get<toml::String>(iter->at("potential"));
        std::unique_ptr<mjolnir::GlobalPotentialBase<traitsT>> pot = 
            read_global_potential<traitsT>(potential, *iter);

        const std::string interaction = 
            toml::get<toml::String>(iter->at("interaction"));
        std::unique_ptr<mjolnir::GlobalInteractionBase<traitsT>> inter = 
            read_global_interaction<traitsT>(interaction, pot, *iter);

        gff.emplace(std::move(inter), std::move(pot));
    }
    return gff;
}

template<typename traitsT>
mjolnir::ForceField<traitsT>
read_force_field(const toml::Table& tab)
{
    std::cerr << "start reading local force field" << std::endl;
    const toml::Array<toml::Table> lff = 
        toml::get<toml::Array<toml::Table>>(tab.at("localforcefield"));
    std::cerr << "end reading local force field" << std::endl;
    const toml::Array<toml::Table> gff = 
        toml::get<toml::Array<toml::Table>>(tab.at("globalforcefield"));
    std::cerr << "end reading global force field" << std::endl;

    mjolnir::ForceField<traitsT> ff(read_local_force_field<traitsT>(lff),
            read_global_force_field<traitsT>(gff));
    return ff;
}

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
    auto pcon   = read_particles<traits>(sim);
    auto rng    = read_random_number_generator<traits>(sim);
    auto integr = read_time_integrator<traits>(sim, rng);
    integr->initialize(pcon); // simulator is responsible for this...?
    const std::size_t total_step = toml::get<toml::Integer>(sim.at("total_step"));
    const std::size_t save_step  = toml::get<toml::Integer>(sim.at("save_step"));

    const std::string filename = toml::get<toml::String>(sim.at("file_name"));
    const std::string trajname = filename + ".xyz";
    const std::string ene_name = filename + ".ene";
    mjolnir::Observer<traits> obs(trajname, ene_name);

    std::cerr << "end reading [simulator] block" << std::endl;

//     auto forcefield = toml::get<toml::Table>(input.at("forcefield"));
    mjolnir::ForceField<traits> ff = read_force_field<traits>(input);

    std::cerr << "end reading [forcefield] block" << std::endl;

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
