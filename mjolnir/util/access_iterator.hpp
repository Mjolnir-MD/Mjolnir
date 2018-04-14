#ifndef MJOLNIR_UTIL_ACCESS_ITERATOR
#define MJOLNIR_UTIL_ACCESS_ITERATOR
#include <type_traits>
#include <iterator>

namespace mjolnir
{

// access member of Iterator::value_type with Accessor
// ex)
// std::vector<std::pair<int, double>> vec = {{1, 3.14}, {2, 6.28}, {3, 9.42}};
// auto accessor = [](const std::pair<int, double>& p) noexcept -> int const& {
//     return p.first;
// };
// for(auto i = make_transform_iterator(vec.cbegin(), accessor),
//          e = make_transform_iterator(vec.cend(), accessor); i!=e; ++i)
// {
//     std::cout << *i << ','; // outputs 1, 2, 3,
// }
template<typename Iterator, typename Accessor>
struct access_iterator
{
    typedef access_iterator<Iterator, Accessor> self_type;
    typedef Iterator base_type;
    typedef Accessor accessor_type;

    typedef std::iterator_traits<base_type> base_traits;
    typedef typename base_traits::difference_type   difference_type;
    typedef typename base_traits::iterator_category iterator_category;

    static constexpr bool is_random_access = std::is_convertible<
        iterator_category, std::random_access_iterator_tag>::value;
    static constexpr bool is_bidirectional = std::is_convertible<
        iterator_category, std::bidirectional_iterator_tag>::value;

    // XXX: this code is c++11, not c++17.
    typedef typename std::result_of<
        accessor_type(typename base_traits::reference)>::type result_type;
    static_assert(std::is_lvalue_reference<result_type>::value,
        "access_iterator requires Accessor to return lvalue reference.");
    typedef result_type reference;

    typedef typename std::remove_const<
        typename std::remove_reference<result_type>::type>::type value_type;
    typedef typename std::add_pointer<
        typename std::remove_reference<result_type>::type>::type pointer;

    access_iterator(const base_type& b, const accessor_type& a)
        : base_(b), accessor_(a)
    {}
    access_iterator(const base_type& b, accessor_type&& a)
        : base_(b), accessor_(std::move(a))
    {}
    ~access_iterator() = default;
    access_iterator(access_iterator const&) = default;
    access_iterator(access_iterator&&)      = default;
    access_iterator& operator=(access_iterator const&) = default;
    access_iterator& operator=(access_iterator&&)      = default;

    reference operator* () const noexcept {return accessor_(*base_);}
    pointer   operator->() const noexcept {return std::addressof(*(*this));}

    self_type& operator++()    noexcept {++base_; return *this;}
    self_type  operator++(int) noexcept
    {const auto tmp(*this); ++base_; return tmp;}

    typename std::enable_if<is_bidirectional, self_type>::type&
    operator--()    noexcept {--base_; return *this;}
    typename std::enable_if<is_bidirectional, self_type>::type
    operator--(int) noexcept {const auto tmp(*this); --base_; return tmp;}

    typename std::enable_if<is_random_access, self_type>::type&
    operator+=(const difference_type d) noexcept {base_ += d; return *this;}
    typename std::enable_if<is_random_access, self_type>::type&
    operator-=(const difference_type d) noexcept {base_ -= d; return *this;}

    typename std::enable_if<is_random_access, difference_type>::type
    operator-(const self_type rhs) noexcept {return this->base_ - rhs.base_;}

    bool operator==(const self_type& rhs) const noexcept
    {return this->base_ == rhs.base_;}
    bool operator!=(const self_type& rhs) const noexcept
    {return this->base_ != rhs.base_;}

    typename std::enable_if<is_random_access, bool>::type
    operator<(const self_type& rhs) const noexcept
    {return this->base_ < rhs.base_;}
    typename std::enable_if<is_random_access, bool>::type
    operator<=(const self_type& rhs) const noexcept
    {return this->base_ <= rhs.base_;}
    typename std::enable_if<is_random_access, bool>::type
    operator>(const self_type& rhs) const noexcept
    {return this->base_ > rhs.base_;}
    typename std::enable_if<is_random_access, bool>::type
    operator>=(const self_type& rhs) const noexcept
    {return this->base_ >= rhs.base_;}

  private:
    base_type     base_;
    accessor_type accessor_;
};
template<typename Iterator, typename Accessor>
constexpr bool access_iterator<Iterator, Accessor>::is_random_access;
template<typename Iterator, typename Accessor>
constexpr bool access_iterator<Iterator, Accessor>::is_bidirectional;

template<typename Iterator, typename Accessor>
inline access_iterator<Iterator, Accessor>
make_access_iterator(Iterator iter, Accessor acc)
{
    return access_iterator<Iterator, Accessor>(iter, acc);
}

} // mjolnir
#endif //MJOLNIR_UTIL_ACCESS_ITERATOR
