+++
title = "External"
weight = 3000
bookCollapseSection = true
+++

# ExternalForceField

`ExternalForceField` contains a set of interactions between particle and an external field.

## [PositionRestraintInteraction]({{<relref "PositionRestraintInteraction.md">}})

It restrains a particle at a fixed point.

It also can restrain a particle at a fixed distance from a certain point.

## [RectangularBoxInteraction]({{<relref "rectangularbox">}})

It wraps the system by a rectangular box. It does not work under the periodic boundary condition.

## [ExternalDistanceInteraction]({{<relref "distance">}})

It depends on the distance between particles and a structure such as a planner surface.

## [AFMFittingInteraction]({{<relref "AFMFittingInteraction.md">}})

It depends on the correlation coefficient between pseudo AFM image and the reference image.

It is capable of fitting a biomolecular model to an experimental AFM image.
