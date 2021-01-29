+++
title = "Integrators"
weight = 400
bookCollapseSection = true
+++

# Integrator

## [G-JFLangevin]({{<relref "G-JFLangevinIntegrator.md">}})

It performs a NVT simulation according to Langevin equation.

It is developed by the following paper.

- [Niels Grønbech-Jensen & Oded Farago, (2013) Mol.Phys. 111:8, 983-991](https://doi.org/10.1080/00268976.2012.760055)

## [BAOABLangevin]({{<relref "BAOABLangevinIntegrator.md">}})

It performs a NVT simulation according to Langevin equation.

It is developed by the following paper.

- [Benedict Leimkuhler and Charles Matthews. Appl. Math. Res. Exp. (2013) 2013:1, pp. 34-56](https://doi.org/10.1093/amrx/abs010)
- [Benedict Leimkuhler and Charles Matthews. J. Chem. Phys. (2013) 138:17, 174102](https://doi.org/10.1063/1.4802990)

## [g-BAOABLangevin]({{<relref "gBAOABLangevinIntegrator.md">}})

It performs a NVT simulation according to Langevin equation.
It is capable of handling constraints.

It is developed by the following paper.

- [Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016) 472:20160138](https://doi.org/10.1098/rspa.2016.0138)

## [UnderdampedLangevin]({{<relref "UnderdampedLangevinIntegrator.md">}})

It performs a NVT simulation according to Langevin equation.

It is introduced in the following papers.

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

It is the same method that is employed in CafeMol.

## [VelocityVerlet]({{<relref "VelocityVerletIntegrator.md">}})

It performs constant energy simulation according to Newtonian equation.
