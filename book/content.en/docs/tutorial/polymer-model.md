+++
title  = "A simple polymer model"
weight = 200
+++

# A simple polymer model

This time, we will run a molecular dynamics simulation of 100 particles connected by springs.

## `[files]` and `[units]`

Same as the previous example, we first define the names and the format of output files.

```toml
[files]
output.prefix = "polymer-model"
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

## `[simulator]`

Unlike the previous example, here, all the particles are connected by springs, so we can omit periodic boundary without worrying that some of the particles will go far away.

```toml
[simulator]
type          = "MolecularDynamics"
precision     = "double"
boundary_type = "Unlimited" # No periodic boundary
seed          = 123456789
delta_t       = 0.01
total_step    = 1_000_000
save_step     =     1_000
integrator.type = "BAOABLangevin"
integrator.gammas = [
    {index =  0, gamma = 1.00},
    {index =  1, gamma = 1.00},
    {index =  2, gamma = 1.00},
    # ...
    {index = 99, gamma = 1.00},
]
```

Most of the parts are the same as the previous example.

To see the complete list of available options, see [Simulator section in the reference.]({{<ref "/docs/reference/simulators">}})

## `[[systems]]`

This time, we don't need `boundary_shape` because we don't have a periodic boundary.

Let's place particles along a straight line.

```toml
[[systems]]
attributes.temperature = 300.0 # K
particles = [
    {mass = 1.0, position = [ 0.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 1.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 2.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 3.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 4.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 5.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 6.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 7.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 8.000, 0.000, 0.000]},
    # 続く...
    {mass = 1.0, position = [99.000, 0.000, 0.000]},
]
```

To see the complete list of available options, see [Systems section in the reference.]({{<ref "/docs/reference/system">}})

## `[[forcefields]]`

Finally, we will define forcefield parameters.

In this example, we want to connect particles by springs. So we define a "local" forcefield.

```toml
[[forcefields]]
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
```

`topology`, later it will be explained in detail, is a label to handle relationships between local and global interactions.
In this case, the name is not important.

Then, apply local potentials between particles.
`BondLength` interaction takes two particles.
`Harmonic` potential takes `v0` and `k`.
`v0` is the value where the potential takes the minimum.
`k` is the value which determines the strength of the potential.

```toml
[[forcefields]]
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    {indices = [ 0, 1], v0 = 1.0, k = 10.0},
    {indices = [ 1, 2], v0 = 1.0, k = 10.0},
    {indices = [ 2, 3], v0 = 1.0, k = 10.0},
    # ...
    {indices = [98,99], v0 = 1.0, k = 10.0},
]
```

To see the complete list of available options, see [ForceFields section in the reference.]({{<ref "/docs/reference/forcefields">}})

## Simulation

Congratulations! All the required information is written in the input file.
By passing the input file to mjolnir, it runs the simulation.

It takes relatively short time compared to the previous example because we have less particles.

```console
$ ./bin/mjolnir polymer-model.toml
reading input file...
-- reading and parsing toml file `polymer-model.toml` ...  successfully parsed.
-- the log file is `./polymer-model.log`
-- mjolnir version v1.22.0-dev (06d5ede6)
-- compiled using /usr/bin/g++-10
-- input_path is ./
-- expanding include files ...
-- done.
-- precision is double
-- Boundary Condition is Unlimited
-- execute on single core
-- energy unit is [kcal/mol]
-- length unit is [angstrom]
-- Integrator is BAOABLangevin.
-- Simulator type is MolecularDynamics.
-- total step is 1000000
-- save step is 1000
-- checkpoint step is 1000
-- reading 0th system ...
-- 100 particles are found.
-- output file prefix is `./polymer-model`
-- output xyz format.
-- LocalForceField (x1) found
-- reading 0th [[forcefields.local]]
-- Bond Length interaction found.
-- -- potential function is Harmonic.
-- -- 99 interactions are found.
-- seed is 123456789
done.
initializing simulator...
-- generating velocity with T = 300...
-- done.
done.
start running simulation
  8.2%|████                                              | 10.5 seconds remaining
```

You will find the following files after the simulation is complete.

```console
$ ls polymer-model*
polymer-model.toml
polymer-model.ene
polymer-model.log
polymer-model_rng.msg
polymer-model_system.msg
polymer-model_position.xyz
polymer-model_velocity.xyz
```

The `.msg` files are for restarting simulation. The format is [MsgPack](https://msgpack.org).

The `.ene` file has the value of energies and other physical quantities in a simple ASCII format, and you can easily visualize it by using gnuplot or other software or a library.

```console
$ head polymer-model.ene
# unit of length : angstrom, unit of energy : kcal/mol
# timestep  BondLength:Harmonic  kinetic_energy attribute:temperature
0                      0.000000      96.956625            300.000000
1000                  56.910239      86.124451            300.000000
2000                  51.497857      85.960034            300.000000
3000                  44.394914      90.433666            300.000000
4000                  50.774066      89.091169            300.000000
5000                  46.749105     104.594276            300.000000
6000                  34.130104      86.045558            300.000000
7000                  39.505721      84.893861            300.000000
```

`polymer-model_position.xyz` contains the trajectory.
You can visualize it by passing the file to a VMD or molecular viewer.

## Conclusion

This tutorial is over.

We passed two particles to each `BondLegnth` interaction.
Those two are not necessarily be contiguous.
You can apply potential on any pair of particles.

`BondAngle` and `DihedralAngle` is also a `LocalInteraction`.
Mjolnir supports those in the same way as `BondLength`.
If you are interested, see [LocalForceField]({{<ref "docs/reference/forcefields/local">}}).
