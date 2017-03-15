#ifndef MEGINGJORD_SIMD_PACKABLE_ARRAY
#define MEGINGJORD_SIMD_PACKABLE_ARRAY

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT use this header alone"
#endif

#include <meginrjord/util/aligned_array.hpp>

namespace megingjord
{
namespace simd
{

/*! @brief aligned array with 0 padding */
template<typename T, std::size_t N, typename simd_traits = MEGINGJORD_DEFAULT_SIMD>
struct packable_array
{
    static_assert(N != 0, "packable_array must have size");
    static_assert(is_packable<T>::value, "packable_array contains only packable type");

    typedef packable_array<T, N, simd_traits> self_type;

    typedef T value_type;
    typedef typename simd_traits::template pack_trait<value_type>::type pack_type;
    typedef typename pack_type::type packed_type;

    constexpr static std::size_t actual_size = N;
    constexpr static std::size_t pack_size   = pack_type::size;
    constexpr static std::size_t align_byte  = pack_type::align_byte;
    constexpr static bool        filled      = (actual_size % pack_size == 0);
    constexpr static std::size_t number_of_packs =
        filled ? (actual_size / pack_size) : (actual_size / pack_size) + 1;
    constexpr static std::size_t filled_size = number_of_packs * pack_size;

    typedef aligned_array<value_type, filled_size, align_byte> container_type;
    typedef typename container_type::pointer                pointer;
    typedef typename container_type::const_pointer          const_pointer;
    typedef typename container_type::reference              reference;
    typedef typename container_type::const_reference        const_reference;
    typedef typename container_type::iterator               iterator;
    typedef typename container_type::const_iterator         const_iterator;
    typedef typename container_type::size_type              size_type;
    typedef typename container_type::difference_type        difference_type;
    typedef typename container_type::reverse_iterator       reverse_iterator;
    typedef typename container_type::const_reverse_iterator const_reverse_iterator;

    container_type values;

    packable_array()  = default;
    ~packable_array() = default;
    packable_array(const packable_array&) = default;
    packable_array(packable_array&&)      = default;
    packable_array& operator=(const packable_array&) = default;
    packable_array& operator=(packable_array&&)      = default;

    packable_array(std::initializer_list<value_type> init);

    constexpr size_type size()     const noexcept {return actual_size;}
    constexpr size_type max_size() const noexcept {return filled_size;}
    constexpr bool      empty()    const noexcept {return false;}
    void fill(const value_type& v){values.fill(v);}

    reference       operator[](const size_type i)       noexcept {return values[i];}
    const_reference operator[](const size_type i) const noexcept {return values[i];}
    reference       at(const size_type i)       {return values.at(i);}
    const_reference at(const size_type i) const {return values.at(i);}

    reference       front()       noexcept {return values.front();}
    const_reference front() const noexcept {return values.front();}
    reference       back()        noexcept {return values.back();}
    const_reference back()  const noexcept {return values.back();}
    pointer         data()        noexcept {return values.data();}
    const_pointer   data()  const noexcept {return values.data();}

    iterator       begin()        noexcept {return values.begin();}
    iterator       end()          noexcept {return values.end();}
    const_iterator begin()  const noexcept {return values.begin();}
    const_iterator end()    const noexcept {return values.end();}
    const_iterator cbegin() const noexcept {return values.cbegin();}
    const_iterator cend()   const noexcept {return values.cend();}

    reverse_iterator       rbegin()        noexcept {return values.rbegin();}
    reverse_iterator       rend()          noexcept {return values.rend();}
    const_reverse_iterator rbegin()  const noexcept {return values.rbegin();}
    const_reverse_iterator rend()    const noexcept {return values.rend();}
    const_reverse_iterator crbegin() const noexcept {return values.crbegin();}
    const_reverse_iterator crend()   const noexcept {return values.crend();}
};

template<typename T, std::size_t N, typename S>
packable_array<T, N, S>::packable_array(std::initializer_list<value_type> init)
{
    if(init.size() > this->max_size()) throw std::out_of_range(
            "packable_array::ctor: size of initializer list exceeds max_size");
    this->fill(0);
    auto v = this->data();
    for(auto iter = init.begin(); iter != init.end(); ++iter)
    {
        *v = std::move(*iter); ++v;
    }
}

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_PACKABLE_ARRAY */
