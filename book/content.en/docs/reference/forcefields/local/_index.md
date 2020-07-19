+++
title  = "Local"
weight = 1000
bookCollapseSection = true
+++

# LocalForceField

`LocalForceField` contains a set of interactions that affects specific particles.

## [BondLengthInteraction]({{<relref "bondlength">}})

It depends on the distance between two particles.

## [BondAngleInteraction]({{<relref "bondangle">}})

It depends on the angle formed by 3 particles.

## [DihedralAngleInteraction]({{<relref "dihedral">}})

It depends on the angle formed by 2 planes, where each plane is defined by 3 particles.

## [ContactInteraction]({{<relref "contact">}})

It depends on the distance between two particles.

Basically it is the same as `BondLength`.
But it consider the case when two particles are too distant and exceeds the cutoff range.

## [DummyInteraction]({{<relref "DummyInteraction.md">}})

It does not do anything. It is used to define topology without applying force.

See [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}) for detail.

## [3SPN2BaseStacking]({{<relref "3SPN2BaseStackingInteraction.md">}})

It is specific to 3SPN2 Coarse-Grained DNA model.
