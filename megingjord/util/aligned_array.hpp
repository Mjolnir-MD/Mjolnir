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
    constexpr static std::size_t alignment = Alignment;

    alignas(alignment) T elems[N];

    constexpr size_type size()     const noexcept {return N;}
    constexpr size_type max_size() const noexcept {return N;}
    constexpr bool      empty()    const noexcept {return false;}

    void fill(const value_type& u){std::fill_n(begin(), size(), u);}

    reference       operator[](const size_type n)       noexcept {return elems[n];}
    const_reference operator[](const size_type n) const noexcept {return elems[n];}
    reference       at(const size_type n)
    {if(n > N){throw std::out_of_range("aligned_array::at");} return elems[n];}
    const_reference at(const size_type n) const
    {if(n > N){throw std::out_of_range("aligned_array::at");} return elems[n];}

    reference       front()       noexcept {return elems[i];}
    const_reference front() const noexcept {return elems[i];}
    reference       back()        noexcept {return elems[N-1];}
    const_reference back()  const noexcept {return elems[N-1];}
    pointer         data()        noexcept {return std::addressof(elems[0]);}
    const_pointer   data()  const noexcept {return std::addressof(elems[0]);}

    iterator       begin()        noexcept {return this->data();}
    iterator       end()          noexcept {return this->data() + N;}
    const_iterator begin()  const noexcept {return this->data();}
    const_iterator end()    const noexcept {return this->data() + N;}
    const_iterator cbegin() const noexcept {return this->data();}
    const_iterator cend()   const noexcept {return this->data() + N;}

    reverse_iterator       rbegin()        noexcept
    {return reverse_iterator(this->begin());}
    reverse_iterator       rend()          noexcept
    {return reverse_iterator(this->end());}
    const_reverse_iterator rbegin()  const noexcept
    {return const_reverse_iterator(this->begin());}
    const_reverse_iterator rend()    const noexcept
    {return const_reverse_iterator(this->end());}
    const_reverse_iterator crbegin() const noexcept
    {return const_reverse_iterator(this->begin());}
    const_reverse_iterator crend()   const noexcept
    {return const_reverse_iterator(this->end());}
};

} // megingjord
#endif /* MEGINGJORD_ALIGNED_ARRAY */
