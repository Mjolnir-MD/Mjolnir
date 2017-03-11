#ifndef MEGINGJORD_SIMD_PACKED_ARRAY
#define MEGINGJORD_SIMD_PACKED_ARRAY

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT use this header alone"
#endif

#include <array>

namespace megingjord
{
namespace simd
{

template<typename T, std::size_t N, typename simd_traits = MEGINGJORD_DEFAULT_SIMD>
struct packed_array
{
    static_assert(N != 0, "packed_array must have size");

    typedef T value_type;
    typedef typename simd_traits::template pack_trait<value_type>::type pack_type;
    typedef typename pack_type::type packed_type;

    constexpr static std::size_t pack_size   = pack_type::size;
    constexpr static std::size_t align_byte  = pack_type::align_byte;
    constexpr static bool        filled      = (N % pack_size == 0);
    constexpr static std::size_t container_size = filled ?
                                         (N / pack_size) : (N / pack_size)+1;
    constexpr static std::size_t filled_size = container_size * pack_size;

    typedef aligned_array<value_type, filled_size, align_byte> aligned_array_type;

    typedef std::array<packed_type, container_size>         container_type;
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

    packed_array()  = default;
    ~packed_array() = default;
    packed_array(const packed_array&) = default;
    packed_array(packed_array&&)      = default;
    packed_array& operator=(const packed_array&) = default;
    packed_array& operator=(packed_array&&)      = default;

    packed_array(const aligned_array_type& xs);
    packed_array& operator=(const aligned_array_type& xs);
    packed_array(const value_type x);

    constexpr size_type size()     const noexcept {return container_size;}
    constexpr size_type max_size() const noexcept {return container_size;}
    constexpr bool      empty()    const noexcept {return false;}

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
packed_array<T, N, S>::packed_array(const aligned_array_type& xs)
{
    const T* ptr = xs.data();
    for(std::size_t i=0; i<container_size; ++i)
    {
        values[i] = load_impl<packed_type>::invoke(ptr);
        ptr += pack_size;
    }
}

template<typename T, std::size_t N, typename S>
packed_array<T, N, S>&
packed_array<T, N, S>::operator=(const aligned_array_type& xs)
{
    const T* ptr = xs.data();
    for(std::size_t i=0; i<container_size; ++i)
    {
        values[i] = load_impl<packed_type>::invoke(ptr);
        ptr += pack_size;
    }
    return *this;
}

template<typename T, std::size_t N, typename S>
packed_array<T, N, S>::packed_array(const value_type x)
{
    for(std::size_t i=0; i<container_size; ++i)
        values[i] = set(x);
}

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_PACKED_ARRAY */
