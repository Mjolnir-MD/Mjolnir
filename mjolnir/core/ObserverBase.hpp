#ifndef MJOLNIR_CORE_OBSERVER_BASE_HPP
#define MJOLNIR_CORE_OBSERVER_BASE_HPP
#include <mjolnir/core/System.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <string>

namespace mjolnir
{

template<typename traitsT>
class ObserverBase
{
  public:
    using traits_type     = traitsT;
    using real_type       = typename traits_type::real_type;
    using coordinate_type = typename traits_type::coordinate_type;
    using system_type     = System<traits_type>;
    using forcefield_type = ForceField<traits_type>;

  public:
    ObserverBase() = default;
    virtual ~ObserverBase() = default;

    virtual void initialize(const std::size_t total_step, const real_type dt,
                            const system_type&, const forcefield_type&) = 0;
    virtual void output    (const std::size_t step,       const real_type dt,
                            const system_type&, const forcefield_type&) = 0;
    virtual void finalize  (const std::size_t total_step, const real_type dt,
                            const system_type&, const forcefield_type&) = 0;

    // for testing purpose.
    virtual std::string const& prefix() const noexcept = 0;
};

namespace detail
{
// it is a helper function to write value as an array of bytes
template<typename T>
void write_as_bytes(std::ostream& os, const T& v) noexcept
{
    using Type = typename std::remove_reference<T>::type;
    os.write(reinterpret_cast<const char*>(std::addressof(v)), sizeof(Type));
    return;
}
} // detail

} // mjolnir
#endif// MJOLNIR_CORE_OBSERVER_BASE_HPP
