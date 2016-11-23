#ifndef MJOLNIR_ZIP_ITERATOR
#define MJOLNIR_ZIP_ITERATOR
#include <tuple>

namespace mjolnir
{

template<typename ... T_args>
class zip_iterator
{
  public:
    using self_type       = zip_iterator<T_args...>;
    using container_type  = std::tuple<T_args...>;
    using container_ptr   = container_type*;
    using container_cptr  = container_type const*;
    using container_ref   = container_type&;
    using container_cref  = container_type const&;

  public:
    zip_iterator(T_args&& ... args): iters_(std::forward(args)...){}
    zip_iterator(const T_args& ... args): iters_(args...){}
    ~zip_iterator() = default;

    self_type& operator++();
    self_type& operator++(int);
    self_type& operator--();
    self_type& operator--(int);

    bool operator==(const self_type& rhs) const;
    bool operator!=(const self_type& rhs) const;

    template<std::size_t I, typename... Types>
    friend typename std::tuple_element<I, std::tuple<Types...>>::type&
    get(zip_iterator<Types...>& t) noexcept;

    template<std::size_t I, typename... Types>
    friend typename std::tuple_element<I, std::tuple<Types...>>::type const&
    get(const zip_iterator<Types...>& t) noexcept;

  private:

    container_type iters_;
};

template<std::size_t I, typename... Types>
typename std::tuple_element<I, std::tuple<Types...>>::type&
get(zip_iterator<Types...>& t) noexcept
{
    return std::get<I>(t.iters_);
}

template<std::size_t I, typename... Types>
typename std::tuple_element<I, std::tuple<Types...>>::type const&
get(const zip_iterator<Types...>& t) noexcept
{
    return std::get<I>(t.iters_);
}

namespace detail
{

template<std::size_t i, typename ... T_args>
struct increment_impl
{
    static std::tuple<T_args...>& advance_front(std::tuple<T_args...>& iters)
    {
        ++(std::get<i-1>(iters));
        return increment_impl<i-1, T_args...>::advance_front(iters);
    }

    static std::tuple<T_args...>& advance_back(std::tuple<T_args...>& iters)
    {
        (std::get<i-1>(iters))++;
        return increment_impl<i-1, T_args...>::advance_back(iters);
    }
};

template<typename ... T_args>
struct increment_impl<0, T_args...>
{
    static std::tuple<T_args...>& advance_front(std::tuple<T_args...>& iters)
    {
        return iters;
    }

    static std::tuple<T_args...>& advance_back(std::tuple<T_args...>& iters)
    {
        return iters;
    }
};

template<std::size_t i, typename ... T_args>
struct decrement_impl
{
    static std::tuple<T_args...>& retrace_front(std::tuple<T_args...>& iters)
    {
        --(std::get<i-1>(iters));
        return decrement_impl<i-1, T_args...>::retrace_front(iters);
    }

    static std::tuple<T_args...>& retrace_back(std::tuple<T_args...>& iters)
    {
        (std::get<i-1>(iters))--;
        return decrement_impl<i-1, T_args...>::retrace_back(iters);
    }
};

template<typename ... T_args>
struct decrement_impl<0, T_args...>
{
    static std::tuple<T_args...>& retrace_front(std::tuple<T_args...>& iters)
    {
        return iters;
    }

    static std::tuple<T_args...>& retrace_back(std::tuple<T_args...>& iters)
    {
        return iters;
    }
};

}// detail

template<typename ... T_args>
inline typename zip_iterator<T_args...>::self_type&
zip_iterator<T_args...>::operator++()
{
    this->iters_ = detail::increment_impl<
        sizeof...(T_args), T_args...>::advance_back(iters_);
    return *this;
}

template<typename ... T_args>
inline typename zip_iterator<T_args...>::self_type&
zip_iterator<T_args...>::operator++(int)
{
    this->iters_ = detail::increment_impl<
        sizeof...(T_args), T_args...>::advance_front(iters_);
    return *this;
}

template<typename ... T_args>
inline typename zip_iterator<T_args...>::self_type&
zip_iterator<T_args...>::operator--()
{
    this->iters_ = detail::decrement_impl<
        sizeof...(T_args), T_args...>::retrace_back(iters_);
    return *this;
}

template<typename ... T_args>
inline typename zip_iterator<T_args...>::self_type&
zip_iterator<T_args...>::operator--(int)
{
    this->iters_ = detail::decrement_impl<
        sizeof...(T_args), T_args...>::retrace_front(iters_);
    return *this;
}

template<typename ... T_args>
inline bool zip_iterator<T_args...>::operator==(const self_type& rhs) const
{
    return this->iters_ == rhs.iters_;
}

template<typename ... T_args>
inline bool zip_iterator<T_args...>::operator!=(const self_type& rhs) const
{
    return this->iters_ != rhs.iters_;
}

} // mjolnir
#endif /* MJOLNIR_ZIP_ITERATOR */
