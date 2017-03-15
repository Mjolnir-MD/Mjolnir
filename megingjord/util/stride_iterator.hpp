#ifndef MEGINGJORD_STRIDE_ITERATOR
#define MEGINGJORD_STRIDE_ITERATOR
#include <iterator>

namespace megingjord
{

template<typename iterT, typename containerT, std::size_t strideN>
struct stride_iterator
{
    typedef iterT iterator_type;
    typedef std::iterator_traits<iterator_type> traits;
    typedef typename traits::difference_type    difference_type;
    typedef typename traits::value_type         value_type;
    typedef typename traits::pointer            pointer;
    typedef typename traits::reference          reference;
    typedef typename traits::iterator_category  iterator_category;
    constexpr static std::size_t stride = strideN;

    static_assert(std::is_same<iterator_category,
                               std::random_access_iterator_tag>::value,
                  "stride_iterator takes only random access iterator");

    iterator_type value;

    stride_iterator() noexcept : value(iterator_type()){}
    stride_iterator(iterator_type iter) noexcept : value(iter){}

    reference operator*()  const noexcept {return *value;}
    pointer   operator->() const noexcept {return value;}
    reference operator[](difference_type n) const noexcept;

    iterator_type raw() const noexcept {return value;}

    stride_iterator& operator++()    noexcept;
    stride_iterator  operator++(int) noexcept;
    stride_iterator& operator--()    noexcept;
    stride_iterator  operator--(int) noexcept;

    stride_iterator& operator+=(std::size_t d) noexcept;
    stride_iterator& operator-=(std::size_t d) noexcept;
    stride_iterator  operator+(std::size_t d) noexcept;
    stride_iterator  operator-(std::size_t d) noexcept;
};

template<typename ptrT, typename ctrT, std::size_t S>
inline typename stride_iterator<ptrT, ctrT, S>::reference
stride_iterator<ptrT, ctrT, S>::operator[](difference_type n) const noexcept
{
    return *(value + n);
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>&
stride_iterator<ptrT, ctrT, S>::operator++() noexcept
{
    value += stride;
    return *this;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>
stride_iterator<ptrT, ctrT, S>::operator++(int) noexcept
{
    stride_iterator tmp(value);
    value += stride;
    return tmp;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>&
stride_iterator<ptrT, ctrT, S>::operator--() noexcept
{
    value -= stride;
    return *this;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>
stride_iterator<ptrT, ctrT, S>::operator--(int) noexcept
{
    stride_iterator tmp(value);
    value -= stride;
    return tmp;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>&
stride_iterator<ptrT, ctrT, S>::operator+=(std::size_t d) noexcept
{
    value += stride * d;
    return *this;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>&
stride_iterator<ptrT, ctrT, S>::operator-=(std::size_t d) noexcept
{
    value -= stride * d;
    return *this;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>
stride_iterator<ptrT, ctrT, S>::operator+(std::size_t d) noexcept
{
    return stride_iterator(value + stride * d);
}

template<typename ptrT, typename ctrT, std::size_t S>
inline stride_iterator<ptrT, ctrT, S>
stride_iterator<ptrT, ctrT, S>::operator-(std::size_t d) noexcept
{
    return stride_iterator(value - stride * d);
}

template<typename ptrT, typename ctrT, std::size_t S>
inline bool
operator==(const stride_iterator<ptrT, ctrT, S>& lhs,
           const stride_iterator<ptrT, ctrT, S>& rhs) noexcept
{
    return lhs.value == rhs.value;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline bool
operator!=(const stride_iterator<ptrT, ctrT, S>& lhs,
           const stride_iterator<ptrT, ctrT, S>& rhs) noexcept
{
    return lhs.value != rhs.value;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline bool
operator<(const stride_iterator<ptrT, ctrT, S>& lhs,
          const stride_iterator<ptrT, ctrT, S>& rhs) noexcept
{
    return lhs.value < rhs.value;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline bool
operator<=(const stride_iterator<ptrT, ctrT, S>& lhs,
           const stride_iterator<ptrT, ctrT, S>& rhs) noexcept
{
    return lhs.value <= rhs.value;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline bool
operator>(const stride_iterator<ptrT, ctrT, S>& lhs,
          const stride_iterator<ptrT, ctrT, S>& rhs) noexcept
{
    return lhs.value > rhs.value;
}

template<typename ptrT, typename ctrT, std::size_t S>
inline bool
operator>=(const stride_iterator<ptrT, ctrT, S>& lhs,
           const stride_iterator<ptrT, ctrT, S>& rhs) noexcept
{
    return lhs.value >= rhs.value;
}

} // megingjord
#endif // MEGINGJORD_STRIDE_ITERATOR
