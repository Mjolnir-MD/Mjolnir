#ifndef AMIME_MAKE_ZIP_ITERATOR
#define AMIME_MAKE_ZIP_ITERATOR
#include <amime/util/zip_iterator.hpp>

namespace amime
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

}// amime
#endif /* AMIME_MAKE_ZIP_ITERATOR */
