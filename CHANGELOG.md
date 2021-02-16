# Change Log

## 1.22.0

### Added

- support OpenMP in G-J&F Langevin
- add WormLikeChainOffsetPotential (@yutakasi634 )

### Changed

- make boundary condition support optional in CMake

### Fixed

- fix several tiny mistakes in the reference documents

## 1.21.0

### Added

- RTree as a spatial partitioning
- G-JF Langevin Integrator

### Fixed

- Fix formula of worm like chain (@yutakasi634 )

## 1.20.0

### Added

- WCA potential
- iSoLF potential
- `(SINGLE|DOUBLE)_PRECISION_SUPPORT=(ON|OFF)` to CMake options

### Changed

- speedup NeighborList construction

### Fixed

- Allow negative offset value in global interaction
- Fix link in documentation (@yutakasi634 )

## 1.19.0

### Added
- support OpenMP + g-BAOAB
- output version and compiler information in log

### Changed
- make boundary condition more customizable

## 1.18.0

### Added

- allow list of indices in offset (@yutakasi634)
- enable to load random number generator state from checkpoint file

### Changed

- update toml11 to v3.6.0

### Fixed

- avoid error in an edge case in PeriodicCellList
- stabilize FLP angle potential calculation in an edge case

## 1.17.0

### Changed

- update toml11 to v3.5.0
- speedup EnergyCalculation simulator

### Fixed

- read boundary information correctly in EnergyCalculation simulator

### Misc

- refactor energy format methods

## 1.16.0

### Added

- documents in English

### Changed

- use Hugo instead of gitbook to host documents
- make MultipleBasin faster

## 1.15.1

### Fixed

- enable to include a file in an array of tables
- avoid infinit recursive file inclusion

## 1.15.0

### Added

- g-BAOABIntegrator (@yutakasi634)
- generalized MultipleBasinModel

### Changed

- check initial condition is in the boundary
- make systems.boundary_shape omittable if there is no boundary

### Fixed

- fix CMakeList for the case when only C++ implementation of OpenMP is found

## 1.14.0

### Added

- RectangularBoxInteraction
  - with `LennardJonesWallPotential` and `ExcludedVolumeWallPotential`

### Changed

- rename `integrator.parameter` to `integrator.gamma` in an input file

## 1.13.1

### Fixed

- read `env` correctly while reading `BAOABLangevin` (@yutakasi634).
- add document about restarting via msgpack

## 1.13.0

### Added

- enable to use offset and env in `[simulator]` table (@yutakasi634)
- check parameter duplication in `[simulator]` table (@yutakasi634)
- enable to include toml files in a input file

The included files will be expanded only once. Recursive inclusion does not work.

```toml
include = "some_table.toml"
[another_table]
include = ["parameter1.toml", "parameter2.toml"]
```

## 1.12.0

### Added

- enable to use index offset in the `[[forcefields]]`
- enable to remove net translation/rotation while time integration

### Changed

- update toml11 to v3.4.0

### Misc

- move interaction and potential codes to forcefield/ directory

## 1.11.0

### Added

- enable to save the last snapshot to msgpack and load it as the initial snapshot

### Changed
- simplify OpenMP implementation
- make dihedral calculation slightly efficient

### Refactoring
- move Topology from System to Forcefield

## 1.10.0

### Added
- Add Inverse power potential (@yutakasi634 )
- Add AFMFlexibleFitting

### Changed

- check energy values are finite when outputting a snapshot
- estimate the remaining time in the progress bar
- upgrade toml11 to v3.3.0

## 1.9.0

### Added

- Add `RepulsiveGoContact` and `AttractiveGoContact`
- Add `WormLikeChainPotential` (@yutakasi634)

### Changed

- Update toml11 to v3.2.1
- Warn if splitted input file has several top-level tables
- Notice `ignore.particles_within` status
- Write system attribute status (e.g. reference `temperature`) to `.ene` file

### Misc

- Add API to re-scale neighbor-list margin
- Move FlexibleLocalPotential from general potential directory to forcefield/FLP

## 1.8.0

### Added

- add EnergyCalculationSimulator

### Changed

- warn if undefined group is found in `[[forcefield.global]]`
- warn if invalid key appears in `[[forcefield.global]]`

### Misc

- enable to clone forcefield
- add initialize() to local potentials (@yutakasi634)
- refactor input funstions

## 1.7.1

### Fixed

- fails to run simulation when non-uniform boundary is applied
- add missing line feed in the error message about dihedral potentials
- fix typo in the name of SI unit (@yutakashi)

### Changed

- speedup initialization of FlexibleLocalDihedral potential

## 1.7.0

### Added

- PDNS interaction
- PWMcos interaction
- DirectionalContactInteraction (@yutakashi634)

### Changed

- update toml11 from v3.0.1 to v3.1.0
- internal architectural changes around neighboring list

## 1.6.3

### Fixed

- fix default periodic cell list
  - Global potentials that are applied to all the particles, such as an excluded volume, are not affected

### Added

- add the `-static-intel` flag to the optimization flags when intel compiler is used
- notice if simulation runs on single core.

## 1.6.2

### Fixed

- bug in openmp implementation of CellLists
  - Do __NOT__ use OpenMP and CellList if a potential is not applied to all the particles
  - Global potentials that are applied to all the particles, such as an excluded volume, are not affected

## 1.6.1

### Fixed

- concurrency bug while generating velocity
  - Do __NOT__ use OpenMP stuff in 1.6.0
- add missing docs about potential combination in Dihedral

### Misc

- refactoring

## 1.6.0

### Added

- add 3SPN2C potential for DNA
- add HardCoreExcludedVolumePotential (@yutakashi634)
- enable to generate initial velocity
- allow boundary.type = "Periodic"

### Changed

- change 3SPN2 input format

### Fixed

- suppress incorrect warnings

## 1.5.2

### Fixed

- cutoff length of gaussian contact
  - Do not use `Contact` interaction with `Gaussian` potential in the earlier releases
- fix some typographical errors in the documentation

### Changed

- change compile flags for test codes

## 1.5.1

### Fixed

- add a workaround for intel compiler linker errors
- add a missing newline in an error message

### Added

- warn if unknown keys appear in a toml file

## 1.5.0

### Added

- support SwitchingForceFieldSimulator
- enable to change cut off ratio of global potentials

### Changed

- optimize efficiency of BondLength<GoContact> and Contact<GoContact>
- add warning messages while reading input file
- upgrade toml11 to v3.0.1

### Misc

- refactoring

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
