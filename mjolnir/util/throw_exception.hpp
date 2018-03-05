#ifndef MJOLNIR_THROW_EXCEPTION_H
#define MJOLNIR_THROW_EXCEPTION_H
#include <utility>
#include <sstream>
#include <string>

namespace mjolnir
{
namespace detail
{

inline void concat_args_to_string_impl(std::ostringstream& oss)
{
    return;
}

template<typename Arg1, typename ...Args>
void concat_args_to_string_impl(std::ostringstream& oss,
        Arg1&& arg1, Args&& ... args)
{
    oss << std::forward<Arg1>(arg1);
    concat_args_to_string_impl(oss, std::forward<Args>(args)...);
    return;
}

template<typename ...Args>
std::string concat_args_to_string(Args&& ... args)
{
    std::ostringstream oss;
    concat_args_to_string_impl(oss, std::forward<Args>(args)...);
    return oss.str();
}

} // detail

template<typename Exception, typename ...Args>
[[noreturn]] void throw_exception(Args&& ... args)
{
    throw Exception(detail::concat_args_to_string(std::forward<Args>(args)...));
}

} // mjolnir
#endif// MJOLNIR_THROW_EXCEPTION_H
