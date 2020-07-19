+++
title = "Global"
weight = 2000
bookCollapseSection = true
+++

# GlobalForceField

`GlobalForceField` constins a set of interactions that affects all the particles.

## [GlobalPair]({{<relref "pair">}})

It depends on the distances between two particles.

## [3SPN2BaseBase](3SPN2BaseBaseInteraction.md)

It is specific to 3SPN2 Coarse-Grained DNA model.

## [ProteinDNANonSpecific](ProteinDNANonSpecificInteraction.md)

It is a Coarse-Grained hydrogen bond model.

## Common parts

There are common fields in `[[forcefields.global]]` table.
One specifies the condition when an interaction is ignored.
Another specifies the spatial partitioning method.

### ignore

There are 3 ways to define ignoring pairs.

First, `particles_within` ignores pairs that are connected within N consecutive topological connection.

Second, `group` ignores inter- or intra-group particle pairs. `group` is defined in [System]({{<relref "/docs/reference/system">}}).

Third, `molecule` ignores inter- or intra-molecule particle pairs. `molecule` is defined as a set of particles that are connected to each other via `bond` topology. That means that, `molecule` is a connected component of a graph when we consider a particle as a node and a `bond` connection as an edge.

To define connections on the topology, see [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}).

- `ignore.particles_within`: Table
  - It ignores pairs of particles that are connected via a certain connection within a certain number.
  - The key is the name of the connection, and the value is the number of concecutive connections.
  - For example, `ignore.particles_within.bond = 3` ignores particles that are connected via 1, 2, or 3 consecutive bond connections.
- `ignore.group`: Table
  - It ignores inter- or intra-group pairs.
  - `intra`: Array of Strings
    - It ignores **intra**-group particle pairs.
  - `inter`: Array of Array of Strings
    - It ignores **inter**-group particle pairs.
  - An example later follows.
- `ignore.molecule`: String
  - It ignores inter-/intra-molecule pairs.
  - `"Nothing"`: It does not ignore anything based on molecule (default).
  - `"Self"`: It ignores intra-molecule pairs.
  - `"Others"`: It ignores inter-molecule pairs.

Let's say we have the following system.

```toml
[[systems]]
particles = [
    {... group = "A"}, # particle 0
    {... group = "A"}, #          1
    {... group = "A"}, #          2
     
    {... group = "B"}, # particle 3
    {... group = "B"}, #          4
    {... group = "B"}, #          5
     
    {... group = "C"}, # particle 6
    {... group = "C"}, #          7
    {... group = "C"}, #          8
]
```

The following interaction

```toml
[[forcefields.global]]
interaction = "Pair"
ignore.group.intra = ["A", "B", "C"]
```

ignores interactions between 0-1, 0-2, 1-2, but does not ignore 0-3, 0-6, 3-6 interactions.

Also, the following interaction

```toml
[[forcefields.global]]
interaction = "Pair"
ignore.group.inter = [["A", "B"], ["A", "C"]]
```

ignores interactions between 0-3, 0-4, 0-5, 3-6, 3-7, 3-8, but does not ignore 0-1, 3-4, 3-6 interactions.

### spatial_partition

It provides spatial partitioning algorithm to construct a neighbor list.

- `type`: String
  - The following spatial partitioning methods are available.
  - `"CellList"`
  - `"VerletList"`
- `margin`: Floating
  - The margin in the neighboring list, relative to the cutoff length.
  - It affects the efficiency, but not the accuracy. The most efficient value depends on a potential to be used.
