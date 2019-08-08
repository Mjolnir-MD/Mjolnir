# Change Log

## 1.4.1

### Changed

- consider energy at cutoff in DebyeHuckel
- add SpatialPartitionBase to reduce compilation time
- replace constant values by constexpr functions

## 1.4.0

### Added

- add 3SPN2 forcefield support
- enable to disable inter- and intra-group interactinos
- enalbe to compile files separately to reduce compilation time

### Misc

- refactoring

## 1.3.2

### Fixed

- `ignore.molecule = "Others"` option
  - Do not use this option in the earlier releases
- suppress warnings

### Added

- Improve tests for ExclusionList

## 1.3.1

### Fixed

- total count of molecules in Topology class
- avoid numerical error for edge cases in FlexibleLocalAngle

### Added

- Improve tests for Topology class

## 1.3.0

### Added

- add BAOABLangevinIntegrator
- add CosinePotential for DihedralAngleInteraction
- enable to define variables in input files

### Changed

- upgrade toml11 to 2.4.0
- disallow to use some of non-periodic potential with DihedralAngle
- make special log message colorful
- check isatty before printing a progress bar

## 1.2.2

### Fixed

- VelocityVerletIntegrator for OpenMP implementation (PR #115)
  - Do not use VelocityVerlet with OpenMP in v1.2.0 and v1.2.1
- output path of log files when files.output.path does not have the last '/'

### Misc

- refactoring

## 1.2.1

### Fixed

- incorrect include file in a test code

## 1.2.0

### Added

- parallelism with OpenMP
- documentation (ja)
- Interaction
  - local/ContactInteraction
  - local/DummyInteraction

### Changed

- upgrade toml11 to 2.3.0
- improve efficiency of Global Interactions

### Misc

- refactoring
  - suppress linter's warnings
  - fix some typographic errors
  - simplify CMakeLists.txt
- improve unit tests
  - add another tests for LocalInteractions

## 1.1.0

### Added

- Interaction
  - external/PositionRestraintInteraction
- Potential
  - external/HarmonicRestraintInteraction
- Observer
  - .trr file format
  - .dcd file format

### Changed

- upgrade Boost to 1.69.0
- upgrade toml11 to 2.2.2

### Misc

- refactoring
  - suppress compiler's warnings
  - unify the coding style
- improve unit tests
  - refactor test codes for input readers
  - add some missing tests

## 1.0.0

### Added

- Simulator
  - MolecularDynamics
  - SimulatedAnnealing
  - SteepestDescent
- Integrator
  - VelocityVerlet
  - UnderdampedLangevin
- Interaction
  - BondLength
  - BondAngle
  - DihedralAngle
  - GlobalPair
  - ExternalDistance
- Potential
  - Harmonic
  - Gaussian
  - GoContact
  - FlexibleLocal
  - LennardJones
  - ExcludedVolume
  - DebyeHuckel
