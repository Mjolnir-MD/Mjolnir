#ifndef MJOLNIR_CORE_LOADER_BASE_HPP
#define MJOLNIR_CORE_LOADER_BASE_HPP
#include <fstream>
#include <string>

namespace mjolnir
{

template<typename traitsT>
class System;

template<typename traitsT>
class LoaderBase
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;

  public:
    LoaderBase() = default;
    virtual ~LoaderBase() {}

    // open files, read header, etc.
    virtual void initialize() = 0;

    virtual std::size_t num_particles() const noexcept = 0;
    virtual std::size_t num_frames()    const noexcept = 0;
    virtual bool        is_eof()        const noexcept = 0;

    // load the next snapshot and write it into the system.
    // If there are no snapshot any more, return false.
    // If the number of particles differs from system, throws runtime_error.
    virtual bool load_next(system_type&) = 0;

    // for testing purpose.
    virtual std::string const& filename() const noexcept = 0;
};

} // mjolnir

#ifdef MJOLNIR_SEPARATE_BUILD
#include <mjolnir/core/SimulatorTraits.hpp>
#include <mjolnir/core/BoundaryCondition.hpp>

namespace mjolnir
{
extern template class LoaderBase<SimulatorTraits<double, UnlimitedBoundary>>;
extern template class LoaderBase<SimulatorTraits<float,  UnlimitedBoundary>>;
extern template class LoaderBase<SimulatorTraits<double, CuboidalPeriodicBoundary>>;
extern template class LoaderBase<SimulatorTraits<float,  CuboidalPeriodicBoundary>>;
} // mjolnir
#endif

#endif//MJOLNIR_CORE_LOADER_BASE_HPP
