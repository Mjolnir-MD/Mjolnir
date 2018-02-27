#ifndef MJOLNIR_ZIP_ITERATOR
#define MJOLNIR_ZIP_ITERATOR
#include <mjolnir/util/index_sequence.hpp>
#include <mjolnir/util/is_all.hpp>
#include <type_traits>
#include <iterator>
#include <tuple>

namespace mjolnir
{
namespace detail
{

template<std::size_t i, typename ... T_args>
struct increment_impl
{
    static inline void invoke(std::tuple<T_args...>& iters)
    {
        ++(std::get<i-1>(iters));
        return increment_impl<i-1, T_args...>::invoke(iters);
    }
};

template<typename ... T_args>
struct increment_impl<0, T_args...>
{
    static inline void invoke(std::tuple<T_args...>& iters) {return;}
};

template<std::size_t i, typename T_diff, typename ... T_args>
struct advance_impl
{
    static inline void invoke(std::tuple<T_args...>& iters, const T_diff& d)
    {
        std::get<i-1>(iters) += d;
        return advance_impl<i-1, T_diff, T_args...>::invoke(iters, d);
    }
};

template<typename T_diff, typename ... T_args>
struct advance_impl<0, T_diff, T_args...>
{
    static inline void invoke(std::tuple<T_args...>& iters, const T_diff& d)
    {return;}
};

template<std::size_t i, typename ... T_args>
struct decrement_impl
{
    static inline void invoke(std::tuple<T_args...>& iters)
    {
        --(std::get<i-1>(iters));
        return decrement_impl<i-1, T_args...>::invoke(iters);
    }
};

template<typename ... T_args>
struct decrement_impl<0, T_args...>
{
    static inline void invoke(std::tuple<T_args...>& iters) {return;}
};

template<std::size_t i, typename T_diff, typename ... T_args>
struct retrace_impl
{
    static inline void invoke(std::tuple<T_args...>& iters, const T_diff& d)
    {
        std::get<i-1>(iters) -= d;
        return retrace_impl<i-1, T_diff, T_args...>::invoke(iters, d);
    }
};

template<typename T_diff, typename ... T_args>
struct retrace_impl<0, T_diff, T_args...>
{
    static inline void invoke(std::tuple<T_args...>& iters, const T_diff& d)
    {return;}
};

template<typename ... Iters, std::size_t ... Idx>
constexpr inline
std::tuple<typename std::iterator_traits<Iters>::reference ...>
pack_references(const std::tuple<Iters...>& iters,
                mjolnir::index_sequence<Idx...>)
{
    return {std::forward_as_tuple(std::get<Idx>(iters).operator*()...)};
}

template<typename ... Iters>
constexpr inline
std::tuple<typename std::iterator_traits<Iters>::reference ...>
pack_references(const std::tuple<Iters...>& is)
{
    return pack_references(is, ::mjolnir::make_index_sequence<sizeof...(Iters)>{});
}

template<typename ... Iters, std::size_t ... Idx>
constexpr inline
std::tuple<typename std::iterator_traits<Iters>::pointer ...>
pack_pointers(const std::tuple<Iters...>& iters,
              mjolnir::index_sequence<Idx...>)
{
    return {std::forward_as_tuple(std::get<Idx>(iters).operator->()...)};
}

template<typename ... Iters>
constexpr inline
std::tuple<typename std::iterator_traits<Iters>::pointer ...>
pack_pointers(const std::tuple<Iters...>& is)
{
    return pack_pointers(is, ::mjolnir::make_index_sequence<sizeof...(Iters)>{});
}
}// detail

template<typename ... T_args>
class zip_iterator
{
  public:
    using self_type      = zip_iterator<T_args...>;
    using container_type = std::tuple<T_args...>;

    using value_type = std::tuple<
        typename std::iterator_traits<T_args>::value_type ...>;
    using pointer    = std::tuple<
        typename std::iterator_traits<T_args>::pointer ...>;
    using reference  = std::tuple<
        typename std::iterator_traits<T_args>::reference ...>;
    using difference_type = typename std::common_type<
        typename std::iterator_traits<T_args>::difference_type ...
        >::type;
    using iterator_category = typename std::common_type<
        typename std::iterator_traits<T_args>::iterator_category ...
        >::type;

  public:
    zip_iterator(const T_args& ... args): iters_(args...){}
    zip_iterator(T_args&& ... args)     : iters_(std::move(args)...){}
    ~zip_iterator() = default;

    zip_iterator(const zip_iterator&) = default;
    zip_iterator(zip_iterator&&)      = default;
    zip_iterator& operator=(const zip_iterator&) = default;
    zip_iterator& operator=(zip_iterator&&)      = default;

    constexpr static std::size_t size() noexcept {return sizeof...(T_args);}

    self_type& operator++() noexcept
    {detail::increment_impl<size(), T_args...>::invoke(iters_); return *this;}
    self_type& operator--() noexcept
    {detail::decrement_impl<size(), T_args...>::invoke(iters_); return *this;}

    self_type& operator++(int) {auto tmp = *this; ++(*this); return tmp;}
    self_type& operator--(int) {auto tmp = *this; --(*this); return tmp;}

    typename std::enable_if<std::is_same<
        iterator_category, std::random_access_iterator_tag
        >::value, self_type>::type&
    operator+=(const difference_type d) noexcept
    {
        detail::advance_impl<size(), difference_type, T_args...
            >::invoke(this->iters_, d);
        return *this;
    }

    typename std::enable_if<std::is_same<
        iterator_category, std::random_access_iterator_tag
        >::value, self_type>::type&
    operator-=(const difference_type d) noexcept
    {
        detail::retrace_impl<size(), difference_type, T_args...
            >::invoke(this->iters_, d);
        return *this;
    }

    reference operator*() const noexcept
    {
        return detail::pack_references(this->iters_);
    }

    pointer operator->() const noexcept
    {
        return detail::pack_pointers(this->iters_);
    }

    template<std::size_t I, typename... Types>
    friend typename std::tuple_element<I, std::tuple<Types...>>::type&
    get(zip_iterator<Types...>& t) noexcept;

    template<std::size_t I, typename... Types>
    friend typename std::tuple_element<I, std::tuple<Types...>>::type const&
    get(const zip_iterator<Types...>& t) noexcept;

    container_type const& base() const noexcept {return iters_;}

  private:

    container_type iters_;
};

template<std::size_t I, typename... Types>
inline typename std::tuple_element<I, std::tuple<Types...>>::type&
get(zip_iterator<Types...>& t) noexcept
{
    return std::get<I>(t.iters_);
}

template<std::size_t I, typename... Types>
inline typename std::tuple_element<I, std::tuple<Types...>>::type const&
get(const zip_iterator<Types...>& t) noexcept
{
    return std::get<I>(t.iters_);
}

template<typename ... Iters>
inline bool operator==(const zip_iterator<Iters...>& lhs,
        const zip_iterator<Iters...>& rhs) noexcept
{
    return lhs.base() == rhs.base();
}

template<typename ... Iters>
inline bool operator!=(const zip_iterator<Iters...>& lhs,
        const zip_iterator<Iters...>& rhs) noexcept
{
    return lhs.base() != rhs.base();
}

template<typename ... Iters>
inline bool operator<(const zip_iterator<Iters...>& lhs,
        const zip_iterator<Iters...>& rhs) noexcept
{
    return lhs.base() < rhs.base();
}

template<typename ... Iters>
inline bool operator<=(const zip_iterator<Iters...>& lhs,
        const zip_iterator<Iters...>& rhs) noexcept
{
    return lhs.base() <= rhs.base();
}

template<typename ... Iters>
inline bool operator>(const zip_iterator<Iters...>& lhs,
        const zip_iterator<Iters...>& rhs) noexcept
{
    return lhs.base() > rhs.base();
}

template<typename ... Iters>
inline bool operator>=(const zip_iterator<Iters...>& lhs,
        const zip_iterator<Iters...>& rhs) noexcept
{
    return lhs.base() >= rhs.base();
}

template<typename ... Iters>
inline zip_iterator<Iters...>
operator+(const zip_iterator<Iters...>& i,
          const typename zip_iterator<Iters...>::difference_type d) noexcept
{
    static_assert(std::is_same<std::random_access_iterator_tag,
        typename zip_iterator<Iters...>::iterator_category>::value,
        "operator+ requires zip_iterator to be random access");
    auto tmp(i); tmp += d;
    return tmp;
}

template<typename ... Iters>
inline zip_iterator<Iters...>
operator+(const typename zip_iterator<Iters...>::difference_type d,
          const zip_iterator<Iters...>& i) noexcept
{
    static_assert(std::is_same<std::random_access_iterator_tag,
        typename zip_iterator<Iters...>::iterator_category>::value,
        "operator+ requires zip_iterator to be random access");
    auto tmp(i); tmp += d;
    return tmp;
}

template<typename ... Iters>
inline typename zip_iterator<Iters...>::difference_type
operator-(const zip_iterator<Iters...>& lhs, const zip_iterator<Iters...>& rhs
          ) noexcept
{
    static_assert(std::is_same<std::random_access_iterator_tag,
        typename zip_iterator<Iters...>::iterator_category>::value,
        "operator- requires zip_iterator to be random access");
    return (get<0>(lhs.base()) - get<0>(rhs.base()));
}

template<typename ... Iters>
inline zip_iterator<Iters...>
operator-(const zip_iterator<Iters...>& i,
          const typename zip_iterator<Iters...>::difference_type d) noexcept
{
    static_assert(std::is_same<std::random_access_iterator_tag,
        typename zip_iterator<Iters...>::iterator_category>::value,
        "operator- requires zip_iterator to be random access");
    auto tmp(i); tmp -= d;
    return tmp;
}

} // mjolnir
#endif /* MJOLNIR_ZIP_ITERATOR */
