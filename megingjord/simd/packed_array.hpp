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

template<typename T, std::size_t N, std::size_t P>
struct packed_array
{
    static_assert(N != 0, "packed_array must have size");
//     static_assert(is_pow_2(P), "pack size must be a power of 2");
    constexpr static std::size_t pack_size = P;
    constexpr static bool filled = N % pack_size == 0;
    constexpr static std::size_t packed_size =
        (filled) ? (N / pack_size) : (N / pack_size)+1;
    constexpr static std::size_t max_size_ = pack_size * packed_size;

    typedef T value_type;
    typedef pack<T, pack_size>                              pack_type;
    typedef typename pack_type::type                        packed_type;
    typedef std::array<packed_type, packed_size>            container_type;
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
    typedef aligned_array<T, N, pack_type::align_byte>      aligned_array_type;
    typedef aligned_array<T, max_size_, pack_type::align_byte> filled_array_type;

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

    constexpr size_type size()     const noexcept {return packed_size;}
    constexpr size_type max_size() const noexcept {return max_size_;}
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

template<typename T, std::size_t N, std::size_t P>
packed_array<T, N, P>::packed_array(const aligned_array_type& xs)
{
    if(filled) //TODO: do this at compile time
    {
        const T* ptr = xs.data();
        for(std::size_t i=0; i<packed_size; ++i)
        {
            values[i] = load_impl<packed_type>::invoke(ptr);
            ptr += pack_size;
        }
    }
    else
    {
        const T* ptr = xs.data();
        for(std::size_t i=0; i<packed_size-1; ++i)
        {
            values[i] = load_impl<packed_type>::invoke(ptr);
            ptr += pack_size;
        }
        std::size_t rest_size = N - (packed_size-1) * pack_size;
        typename pack_type::array_type ar;
        for(std::size_t i=0; i<pack_size; ++i)
        {
            if(i < rest_size)
                ar[i] = xs[(packed_size-1) * pack_size + i];
            else
                ar[i] = 0;
        }
        values[packed_size-1] = load(ar);
    }
}

template<typename T, std::size_t N, std::size_t P>
packed_array<T, N, P>&
packed_array<T, N, P>::operator=(const aligned_array_type& xs)
{
    if(filled) //TODO: do this at compile time
    {
        const T* ptr = xs.data();
        for(std::size_t i=0; i<packed_size; ++i)
        {
            values[i] = load_impl<packed_type>::invoke(ptr);
            ptr += pack_size;
        }
    }
    else
    {
        const T* ptr = xs.data();
        for(std::size_t i=0; i<packed_size-1; ++i)
        {
            values[i] = load_impl<packed_type>::invoke(ptr);
            ptr += pack_size;
        }
        std::size_t rest_size = N - (packed_size-1) * pack_size;
        typename pack_type::array_type ar;
        for(std::size_t i=0; i<pack_size; ++i)
        {
            if(i < rest_size)
                ar[i] = xs[(packed_size-1) * pack_size + i];
            else
                ar[i] = 0;
        }
        values[packed_size-1] = load(ar);
    }
}

template<typename T, std::size_t N, std::size_t P>
packed_array<T, N, P>::packed_array(const value_type x)
{
    if(filled) //TODO: do this at compile time
    {
        for(std::size_t i=0; i<packed_size; ++i)
            values[i] = set(x);
    }
    else
    {
        for(std::size_t i=0; i<packed_size-1; ++i)
            values[i] = set(x);

        std::size_t rest_size = N - (packed_size-1) * pack_size;
        typename pack_type::array_type ar;
        for(std::size_t i=0; i<pack_size; ++i)
        {
            if(i < rest_size)
                ar[i] = x;
            else
                ar[i] = 0;
        }
        values[packed_size-1] = load(ar);
    }
}

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_PACKED_ARRAY */
