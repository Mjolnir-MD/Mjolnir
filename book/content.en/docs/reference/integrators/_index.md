+++
title = "Integrators"
weight = 400
bookCollapseSection = true
+++

# Integrator

## [BAOABLangevin]({{<relref "BAOABLangevinIntegrator.md">}})

It performs a NVT simulation according to Langevin equation.

It is developed by the following paper.

- Leimkuhler B, Matthews C. Appl. Math. Res. Exp. (2013)
- Leimkuhler B, Matthews C. J. Chem. Phys. (2013)

## [g-BAOABLangevin]({{<relref "gBAOABLangevinIntegrator.md">}})

It performs a NVT simulation according to Langevin equation.
It is capable of handling constraints.

It is developed by the following paper.

- Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016)

## [UnderdampedLangevin]({{<relref "UnderdampedLangevinIntegrator.md">}})

It performs a NVT simulation according to Langevin equation.

It is introduced in the following papers.

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

It is the same method that is employed in CafeMol.

## [VelocityVerlet]({{<relref "VelocityVerletIntegrator.md">}})

It performs constant energy simulation according to Newtonian equation.
