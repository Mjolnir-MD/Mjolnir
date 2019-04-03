#ifndef MJOLNIR_UTIL_ALIGNED_STORAGE_HPP
#define MJOLNIR_UTIL_ALIGNED_STORAGE_HPP
#include <type_traits>
#include <array>
#include <cstdint>

namespace mjolnir
{

template<typename T, std::size_t Align = 64>
struct alignas(Align) aligned_storage
{
  public:
    using type = T;
    static_assert(Align != 0, "");

    static constexpr std::size_t alignment_size = Align;
    static constexpr std::size_t value_size     = sizeof(T);
    static constexpr std::size_t padded_size    =
        alignment_size * (value_size / alignment_size + 1);
    static constexpr std::size_t padding_size =
        alignment_size - (value_size % alignment_size);

    static_assert((padding_size + value_size) % padded_size    == 0, "");
    static_assert((padding_size + value_size) % alignment_size == 0, "");

  public:
    aligned_storage() = default;
    ~aligned_storage() = default;
    aligned_storage(aligned_storage const&) = default;
    aligned_storage(aligned_storage &&)     = default;
    aligned_storage& operator=(aligned_storage const&) = default;
    aligned_storage& operator=(aligned_storage &&)     = default;

    aligned_storage(const T& v)
        noexcept(std::is_nothrow_copy_constructible<T>::value)
        : value(v)
    {}
    aligned_storage(T&& v)
        noexcept(std::is_nothrow_move_constructible<T>::value)
        : value(v)
    {}

    T value;

  private:
    std::array<std::uint8_t, padding_size> padding_;
};

static_assert(alignof(aligned_storage<std::int32_t, 64>) == 64, "");
static_assert(sizeof (aligned_storage<std::int32_t, 64>) == 64, "");

} // mjolnir
#endif// MJOLNIR_UTIL_ALIGNED_STORAGE_HPP
