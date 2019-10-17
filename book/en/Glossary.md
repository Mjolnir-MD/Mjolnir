# Glossary

## Simulator

Simulator implements a simulation protocol.

`MolecularDynamicsSimulator` just performs time integration using given integrator, forcefields and systems.

`SimulatedAnnealingSimulator` also performs time integration and changes temperature of the system in the specified manner.

`SwitchingForceFieldSimulator` manages multiple forcefields and swaps the active forcefield at the specified timestep.

The actual force and time integral calculations are performed by another component managed by the Simulator.

## System

A group of particles to be simulated.

In addition to particles, it also manages boundary conditions, system parameters (such as temperature), and Topology.

## Topology

This is a (mathematical) graph showing which particles are connected by which interaction.

It is used in a ForceField to check particles are already connected.

## Integrator

Integrator performs time integration.

## ForceField

There are three types: local, global, and external.

### LocalForceField

It is an interaction between specific particles.

Typical examples are bond length potential, bond angle potential, and dihedral potential.

### GlobalForceField

It is an interaction between all the combinations of particles.

Typical examples are intermolecular force potential and electrostatic potential.

### ExternalForceField

It is an interaction between particles and the external field.

Only this can cause translation and rotation of the entire system.

Typical examples are anchors to specific points in the space and boxes that cover the entire system.

### Interaction

ForceField consists of Interaction and Potential.

Interaction manages on which particle the potential will be applied.

When there is a force field like

{% math%}
U(\theta) = k(\theta-\theta_0)^2,
{% endmath%}

the force becomes the following.

{% math%}
\require{color}
\nabla_i U(\theta) = \color[RGB]{255,0,0} (\nabla_i \theta) \color[RGB]{0,0,0} \frac{dU}{d\theta}
{% endmath%}

It implements the first term that is colored red.

### Potential

ForceField consists of Interaction and Potential.

Potential manages the parameters related to the potential function and calculates the potential function and the derivative.

When there is a force field like

{% math%}
U(\theta) = k (\theta - \theta_0)^2,
{% endmath%}

the force becomes the following.

{% math%}
\require{color}
\nabla_i U(\theta) = \color[RGB]{255,0,0} (\nabla_i \theta) \color[RGB]{0,0,0} \frac{dU}{d\theta}
{% endmath%}

It implements the second term that is colored red.

## Observer

Observer outputs energy, particle position, and velocity to a file.

## Traits

A tag type to dispatch the implementations.

For example, codes needed only for OpenMP are specialized with `OpenMPSimulatorTraits`.
