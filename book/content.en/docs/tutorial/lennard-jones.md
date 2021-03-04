+++
title  = "The Lennard-Jones fluid"
weight = 100
+++

# The Lennard-Jones fluid

As a simple example, we will run a molecular dynamics simulation of particles interacting via Lennard-Jones potential.
We want to make the system small so that the simulation does not takes long time, we put only 512 particles.
Also, our goal is not a simulation of a specific atom but just to run a simulation, so we keep parameters simple.

We need the following information to run a simulation.

1. Simulation method
2. Initial configuration
3. Forcefield parameters

We will write those information and save it in an input file named `lennard-jones.toml`.

## `[files]` and `[units]`

First of all, we will define the output file name and format.
To make it simple, we use [xyz](https://en.wikipedia.org/wiki/XYZ_file_format) format.

```toml
[files]
output.prefix = "lennard-jones"
output.path   = "./"
output.format = "xyz"
```

To see the complete list of available options, see [files section in the reference.]({{<ref "/docs/reference/files">}})

Also, we need to define unit system to be used.

```toml
[units]
length = "angstrom"
energy = "kcal/mol"
```

To see the complete list of available options, see [units section in the reference.]({{<ref "/docs/reference/units">}})

## `[simulator]` and other general properties

We have defined general settings, so next we will define the simulation method.

We will perform normal molecular dynamics simulation, so the type of simulator is `MolecularDynamics`.
We will apply the periodic boundary conditions so that the simulated particles are enclosed in a box.
And, while the simulation, we use double precision floating point.
We also need to set the seed of random number generator.

```toml
[simulator]
type          = "MolecularDynamics"
precision     = "double"
boundary_type = "Periodic"
seed          = 123456789
```

The integrator is also needed.
Here, we use BAOAB-type Langevin integrator that is introduced relatively recently.
We can set the friction coefficient one by one, but for now, we just set 1 for all the particles.

```toml
[simulator]
type          = "MolecularDynamics"
precision     = "double"
boundary_type = "Periodic"
seed          = 123456789
delta_t       = 0.01
total_step    = 100_000
save_step     =     100
integrator.type = "BAOABLangevin"
integrator.gammas = [
    {index =   0, gamma = 1.00},
    {index =   1, gamma = 1.00},
    {index =   2, gamma = 1.00},
    # ...
    {index = 511, gamma = 1.00},
]
```

It is a hard work. But keyboard macro supported by the editor or a simple script can help you to generate it.

To see the complete list of available options, see [Simulator section in the reference.]({{<ref "/docs/reference/simulators">}})

## `[[systems]]`

To run a simulation, we need a reasonable initial configuration.
Initial configuration is defined in `[[systems]]`.

First, we will define a general parameter.
The temperature of the system and the shape of periodic boundary.

```toml
[[systems]]
attributes.temperature = 300.0 # K
boundary_shape.lower  = [ 0.0,  0.0,  0.0]
boundary_shape.upper  = [16.0, 16.0, 16.0]
```

To avoid collisions, we put particles on grid points separating them by 2 angstroms, assuming they have radius of 1 angstrom..

```toml
[[systems]]
attributes.temperature = 300.0 # K
boundary_shape.lower  = [ 0.0,  0.0,  0.0]
boundary_shape.upper  = [16.0, 16.0, 16.0]
particles = [
    {mass = 1.0, position = [ 1.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 3.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 5.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 7.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 9.000, 1.000, 1.000]},
    {mass = 1.0, position = [11.000, 1.000, 1.000]},
    {mass = 1.0, position = [13.000, 1.000, 1.000]},
    {mass = 1.0, position = [15.000, 1.000, 1.000]},

    {mass = 1.0, position = [ 1.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 3.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 5.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 7.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 9.000, 3.000, 1.000]},
    {mass = 1.0, position = [11.000, 3.000, 1.000]},
    {mass = 1.0, position = [13.000, 3.000, 1.000]},
    {mass = 1.0, position = [15.000, 3.000, 1.000]},
    # 続く...
]
```

To see the complete list of available options, see [Systems section in the reference.]({{<ref "/docs/reference/system">}})

## `[[forcefields]]`

Finally, we will define forcefield parameters.

Lennard-Jones interaction is applied to all the pair of particles, so it is non-local interaction.
To distinguish from external force field, we call it as a global interaction.
So first define a `[[forcefields.global]]` under the `[[forcefields]]` table and specify the interaction type and the potential function.

```toml
[[forcefields]]
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJones"
```

Next, we will specify the spatial partitioning method.
In this example, our system could be too small to construct an efficient cell list.
So we use raw [Verlet list](https://en.wikipedia.org/wiki/Verlet_list).

Verlet list takes a margin to reduce the computational cost without changing the simulation result.
By taking particles that is a bit distant from cutoff region account, the list can be reused until a particle that is not registered in a list comes into the cutoff region.

The `margin` is defined relative to the cutoff length.
To tune it, we need to run several short simulations with varying margin.
But since this does not change the trajectory, any value is basically okay (if the simulation doesn't take too long time).

```toml
[[forcefields]]
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJones"
spatial_partition.type = "VerletList"
spatial_partition.margin = 0.25
```

And of course, we need a forcefield parameter for each particle.

```toml
[[forcefields]]
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJones"
spatial_partition.type = "VerletList"
spatial_partition.margin = 0.25
parameters = [
    {index =   0, sigma = 2.0, epsilon = 1.0},
    {index =   1, sigma = 2.0, epsilon = 1.0},
    {index =   2, sigma = 2.0, epsilon = 1.0},
    # ...
    {index = 511, sigma = 2.0, epsilon = 1.0},
]
```

The actual parameter for each pair or particles is calculated according to the Lorentz-Berthelot combining rule.

To see the complete list of available options, see [ForceFields section in the reference.]({{<ref "/docs/reference/forcefields">}})

## Simulation

Congratulations! All the required information is written in the input file.
By passing the input file to mjolnir, it runs the simulation.
It takes a few minutes, depending on the machine you are using.

```console
$ ./bin/mjolnir lennard-jones.toml
reading input file...
-- reading and parsing toml file `lennard-jones.toml` ...  successfully parsed.
-- the log file is `./lennard-jones.log`
-- mjolnir version v1.22.0
-- compiled using /usr/bin/g++-10
-- input_path is ./
-- expanding include files ...
-- done.
-- precision is double
-- Boundary Condition is Periodic. The shape is cuboid.
-- execute on single core
-- energy unit is [kcal/mol]
-- length unit is [angstrom]
-- Integrator is BAOABLangevin.
-- Simulator type is MolecularDynamics.
-- total step is 100000
-- save step is 100
-- checkpoint step is 100
-- reading 0th system ...
-- 512 particles are found.
-- output file prefix is `./lennard-jones`
-- output xyz format.
-- GlobalForceField (x1) found
-- reading 0th [[forcefields.global]]
-- Pair interaction found.
-- -- potential function is Lennard-Jones.
-- -- Spatial Partition is VerletList with relative margin = 0.25
-- No `ignore.group` is provided. All the groups are taken into account.
-- No `ignore.molecule` is provided. All the molecules are taken into account.
-- seed is 123456789
done.
initializing simulator...
-- generating velocity with T = 300...
-- done.
done.
start running simulation
 14.1%|███████                                           | 56.0 seconds remaining
```

You will find the following files after the simulation is complete.

```console
$ ls lennard-jones*
lennard-jones.toml
lennard-jones.log
lennard-jones.ene
lennard-jones_position.xyz
lennard-jones_velocity.xyz
lennard-jones_system.msg
lennard-jones_rng.msg
```

The `.msg` files are for restarting simulation. The format is [MsgPack](https://msgpack.org).

The `.ene` file has the value of energies and other physical quantities in a simple ASCII format, and you can easily visualize it by using gnuplot or other software or a library.

```console
$ head lennard-jones.ene
# unit of length : angstrom, unit of energy : kcal/mol
# timestep  GlobalPairLennardJones  kinetic_energy attribute:temperature
0                     -1704.786330     453.454993            300.000000
100                   -2579.996798     812.198804            300.000000
200                   -2787.935537     554.446648            300.000000
300                   -2909.224864     473.619251            300.000000
400                   -2913.189150     453.464065            300.000000
500                   -2964.234777     463.400579            300.000000
600                   -2988.270454     462.704726            300.000000
700                   -2960.111833     458.723132            300.000000
```

`lennard-jones_position.xyz` contains the trajectory.
You can visualize it by passing the file to a VMD or molecular viewer.

## Conclusion

This tutorial is over.

You might feel the input file is too explicit.
But, it allows you to easily change the forcefield parameter one by one.
You can increase the size of one particular particle, or mix different-sized particles.
If you are interested, have fun with it.
