# Change Log

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
