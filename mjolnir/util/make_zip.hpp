#ifndef MJOLNIR_MAKE_ZIP_ITERATOR
#define MJOLNIR_MAKE_ZIP_ITERATOR
#include "zip_iterator.hpp"
#include <iterator>

namespace mjolnir
{

template<typename ... T_args>
inline zip_iterator<T_args...> make_zip(const T_args& ... args)
{
    return zip_iterator<T_args...>(args...);
}

template<typename ... T_args>
inline zip_iterator<T_args...> make_zip_cbegin(const T_args& ... args)
{
    return zip_iterator<T_args...>(std::begin(args)...);
}

template<typename ... T_args>
inline zip_iterator<T_args...> make_zip_begin(T_args& ... args)
{
    return zip_iterator<T_args...>(std::begin(args)...);
}

template<typename ... T_args>
inline zip_iterator<T_args...> make_zip_cend(const T_args& ... args)
{
    return zip_iterator<T_args...>(std::end(args)...);
}

template<typename ... T_args>
inline zip_iterator<T_args...> make_zip_end(T_args& ... args)
{
    return zip_iterator<T_args...>(std::end(args)...);
}

}// mjolnir
#endif /* MJOLNIR_MAKE_ZIP_ITERATOR */
