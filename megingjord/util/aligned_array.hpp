#ifndef MEGINGJORD_ALIGNED_ARRAY
#define MEGINGJORD_ALIGNED_ARRAY
#include <stdexcept>
#include <algorithm>
#include <cstddef>

namespace megingjord
{

template<typename T, std::size_t N, std::size_t N_align = 32>
struct aligned_array
{
    static_assert(N != 0, "aligned_array must have size");

    typedef T value_type;
    typedef value_type*       pointer;
    typedef value_type const* const_pointer;
    typedef value_type&       reference;
    typedef value_type const& const_reference;
    typedef value_type*       iterator;
    typedef value_type const* const_iterator;
    typedef std::size_t       size_type;
    typedef std::ptrdiff_t    difference_type;
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    constexpr static std::size_t alignment = N_align;

    alignas(alignment) T elems[N];

    constexpr size_type size()     const noexcept {return N;}
    constexpr size_type max_size() const noexcept {return N;}
    constexpr bool      empty()    const noexcept {return false;}

    void fill(const value_type& u){std::fill_n(begin(),size(),u);}

    reference       operator[](const size_type n)       noexcept;
    const_reference operator[](const size_type n) const noexcept;
    reference       at(const size_type n);
    const_reference at(const size_type n) const;

    reference       front()       noexcept;
    const_reference front() const noexcept;
    reference       back()        noexcept;
    const_reference back()  const noexcept;
    pointer         data()        noexcept;
    const_pointer   data()  const noexcept;

    iterator       begin()        noexcept;
    iterator       end()          noexcept;
    const_iterator begin()  const noexcept;
    const_iterator end()    const noexcept;
    const_iterator cbegin() const noexcept;
    const_iterator cend()   const noexcept;

    reverse_iterator       rbegin()        noexcept;
    reverse_iterator       rend()          noexcept;
    const_reverse_iterator rbegin()  const noexcept;
    const_reverse_iterator rend()    const noexcept;
    const_reverse_iterator crbegin() const noexcept;
    const_reverse_iterator crend()   const noexcept;
};

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::reference
aligned_array<T, N, Na>::operator[](const size_type n) noexcept
{
    return elems[n];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reference
aligned_array<T, N, Na>::operator[](const size_type n) const noexcept
{
    return elems[n];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::reference
aligned_array<T, N, Na>::at(const size_type n)
{
    if(n >= N) throw std::out_of_range("aligned_array::at: n >= N");
    return elems[n];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reference
aligned_array<T, N, Na>::at(const size_type n) const
{
    if(n >= N) throw std::out_of_range("aligned_array::at: n >= N");
    return elems[n];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::reference
aligned_array<T, N, Na>::front() noexcept
{
    return elems[0];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reference
aligned_array<T, N, Na>::front() const noexcept
{
    return elems[0];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::reference
aligned_array<T, N, Na>::back() noexcept
{
    return elems[N-1];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reference
aligned_array<T, N, Na>::back() const noexcept
{
    return elems[N-1];
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::pointer
aligned_array<T, N, Na>::data() noexcept
{
    return elems;
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_pointer
aligned_array<T, N, Na>::data() const noexcept
{
    return elems;
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::iterator
aligned_array<T, N, Na>::begin() noexcept
{
    return iterator(data());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::iterator
aligned_array<T, N, Na>::end() noexcept
{
    return iterator(data() + N);
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_iterator
aligned_array<T, N, Na>::begin() const noexcept
{
    return const_iterator(data());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_iterator
aligned_array<T, N, Na>::end() const noexcept
{
    return const_iterator(data() + N);
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_iterator
aligned_array<T, N, Na>::cbegin() const noexcept
{
    return const_iterator(data());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_iterator
aligned_array<T, N, Na>::cend() const noexcept
{
    return const_iterator(data() + N);
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::reverse_iterator
aligned_array<T, N, Na>::rbegin() noexcept
{
    return reverse_iterator(end());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::reverse_iterator
aligned_array<T, N, Na>::rend() noexcept
{
    return reverse_iterator(begin());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reverse_iterator
aligned_array<T, N, Na>::rbegin() const noexcept
{
    return const_reverse_iterator(end());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reverse_iterator
aligned_array<T, N, Na>::rend() const noexcept
{
    return const_reverse_iterator(begin());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reverse_iterator
aligned_array<T, N, Na>::crbegin() const noexcept
{
    return const_reverse_iterator(end());
}

template<typename T, std::size_t N, std::size_t Na>
inline typename aligned_array<T, N, Na>::const_reverse_iterator
aligned_array<T, N, Na>::crend() const noexcept
{
    return const_reverse_iterator(begin());
}

} // megingjord
#endif /* MEGINGJORD_ALIGNED_ARRAY */
