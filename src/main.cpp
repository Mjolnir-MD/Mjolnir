#include <mjolnir/core/Simulator.hpp>
#include <mjolnir/core/Observer.hpp>
#include <mjolnir/core/VelocityVerlet.hpp>
#include <mjolnir/core/GlobalDistanceInteraction.hpp>
#include <mjolnir/core/LennardJonesPotential.hpp>
#include <mjolnir/core/NaivePairCalculation.hpp>
#include <mjolnir/core/VerletList.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>
#include <random>

typedef mjolnir::DefaultTraits traits;

int main()
{
    // parameters
    const std::size_t num_particle = 1000;
    const double radius = 3.0;
    const double dt = 0.01;

    std::mt19937 mt(134573984);
    std::normal_distribution<double> dist(0.0, 0.001);

    // create particles
    mjolnir::ParticleContainer<traits> pcon(num_particle);
    for(std::size_t i=0; i<num_particle; ++i)
    {
        const std::size_t xi = i % 10;
        const std::size_t yi = ((i-xi) / 10) % 10;
        const std::size_t zi = ((i-xi - yi) / 100) % 10;
        const double vx = dist(mt);
        const double vy = dist(mt);
        const double vz = dist(mt);
        const mjolnir::Vector<double, 3> pos(
                std::pow(2, 1./6.) * radius * xi,
                std::pow(2, 1./6.) * radius * yi,
                std::pow(2, 1./6.) * radius * zi);
        const mjolnir::Vector<double, 3> vel(vx, vy, vz);
        const mjolnir::Vector<double, 3> f(0., 0., 0.);

        pcon[i] = mjolnir::make_particle(1.0, pos, vel, f);
    }

    // create empty local forcefield
    mjolnir::LocalForceField<traits> lff;

    // craete global forcefield
    mjolnir::GlobalForceField<traits> gff;

    // add Lennard-Jones to GlobalForceField
    auto lj = mjolnir::make_unique<mjolnir::LennardJonesPotential<traits>>();
    for(std::size_t i=0; i<num_particle; ++i)
        lj->emplace(std::make_pair(radius, 1.0));

    // prepair verlet list
    traits::real_type cutoff = lj->max_cutoff_length();
    auto space = mjolnir::make_unique<mjolnir::VerletList<traits>>(
            cutoff, cutoff * 1.25, dt);
    space->make(pcon);

    auto gint = mjolnir::make_unique<
        mjolnir::GlobalDistanceInteraction<traits>>(std::move(space));

    gff.emplace(std::move(gint), std::move(lj));

    // move local & global forcefield into "forcefield"
    mjolnir::ForceField<traits> ff(std::move(lff), std::move(gff));

    // create velocity verlet time integration method
    auto integr = mjolnir::make_unique<mjolnir::VelocityVerlet<traits>>(
            dt, num_particle);
    integr->initialize(pcon); // simulator is responsible for this...?

    // create simulator with particle_container, forcefield, time_integrator
    mjolnir::Simulator<traits> sim(
            std::move(pcon), std::move(ff), std::move(integr));

    // create observer
    mjolnir::Observer<traits> obs("test.xyz", "test.ene");

    // run md
    for(std::size_t i=0; i<10000; ++i)
    {
        if(i%100 == 0)
        {
            obs.output_coordinate(sim);
            obs.output_energy(sim);
        }
        sim.step();
    }

    // output last state
    obs.output_coordinate(sim);
    obs.output_energy(sim);

    return 0;
}
