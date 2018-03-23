#ifndef MEGINGJORD_STRIDE_ITERATOR
#define MEGINGJORD_STRIDE_ITERATOR
#include <iterator>

namespace megingjord
{

template<typename Iterator, std::size_t Stride>
struct stride_iterator
{
    typedef Iterator iterator_type;
    static constexpr std::size_t stride = Stride;

    typedef std::iterator_traits<iterator_type>    traits;
    typedef typename traits::difference_type       difference_type;
    typedef typename traits::value_type            value_type;
    typedef typename traits::pointer               pointer;
    typedef typename traits::reference             reference;
    typedef typename traits::iterator_category     iterator_category;
    typedef stride_iterator<iterator_type, stride> self_type;

    static_assert(std::is_same<iterator_category,
                               std::random_access_iterator_tag>::value,
                  "stride_iterator takes only random access iterator");

    stride_iterator() = default;
    ~stride_iterator() = default;
    stride_iterator(stride_iterator const&) = default;
    stride_iterator(stride_iterator&&)      = default;
    stride_iterator& operator=(stride_iterator const&) = default;
    stride_iterator& operator=(stride_iterator&&)      = default;

    stride_iterator(const iterator_type iter) noexcept : base_(iter){}

    reference operator* () const noexcept {return base_.operator*();}
    pointer   operator->() const noexcept {return base_.operator->();}
    reference operator[](difference_type n) const noexcept {return base_[n];}

    self_type& operator++()    noexcept {this->base_ += stride; return *this;}
    self_type& operator--()    noexcept {this->base_ -= stride; return *this;}
    self_type  operator++(int) noexcept
    {const auto tmp = *this; ++(*this); return tmp;}
    self_type  operator--(int) noexcept
    {const auto tmp = *this; --(*this); return tmp;}

    stride_iterator& operator+=(const difference_type d) noexcept
    {this->base_ += d * stride; return *this;}
    stride_iterator& operator-=(const difference_type d) noexcept
    {this->base_ -= d * stride; return *this;}

    stride_iterator operator+(const difference_type d) const noexcept;
    {auto tmp = *this; tmp += d; return tmp;}
    stride_iterator operator-(const difference_type d) const noexcept;
    {auto tmp = *this; tmp -= d; return tmp;}

    differnce_type operator-(const self_type rhs) const noexcept;
    {return *this - rhs;}

    iterator_type const& base() const noexcept {return base_;}

  private:

    iterator_type base_;
};

template<typename I, std::size_t S>
constexpr std::size_t stride_iterator<I, S>::stride;

template<typename I, std::size_t S>
inline bool operator==(const stride_iterator<I, S>& lhs,
                       const stride_iterator<I, S>& rhs) noexcept
{
    return lhs.base() == rhs.base();
}

template<typename I, std::size_t S>
inline bool operator!=(const stride_iterator<I, S>& lhs,
                       const stride_iterator<I, S>& rhs) noexcept
{
    return lhs.base() != rhs.base();
}

template<typename I, std::size_t S>
inline bool operator<(const stride_iterator<I, S>& lhs,
                      const stride_iterator<I, S>& rhs) noexcept
{
    return lhs.base() < rhs.base();
}

template<typename I, std::size_t S>
inline bool operator<=(const stride_iterator<I, S>& lhs,
                       const stride_iterator<I, S>& rhs) noexcept
{
    return lhs.base() <= rhs.base();
}

template<typename I, std::size_t S>
inline bool operator>(const stride_iterator<I, S>& lhs,
                      const stride_iterator<I, S>& rhs) noexcept
{
    return lhs.base() > rhs.base();
}

template<typename I, std::size_t S>
inline bool operator>=(const stride_iterator<I, S>& lhs,
                       const stride_iterator<I, S>& rhs) noexcept
{
    return lhs.base() >= rhs.base();
}

} // megingjord
#endif // MEGINGJORD_STRIDE_ITERATOR
