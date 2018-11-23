#ifndef MJOLNIR_UTIL_LOGGER
#define MJOLNIR_UTIL_LOGGER
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <array>
#include <vector>
#include <map>
#include <unordered_map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>

namespace mjolnir
{

/* these macros are defined.                                                 *
 * - MJOLNIR_LOG_DEBUG  -- for debugging. it makes simulation too slow.      *
 * - MJOLNIR_LOG_INFO   -- output in-detail info unless it slows down.       *
 * - MJOLNIR_LOG_NOTICE -- write current progress and status to console.     *
 * - MJOLNIR_LOG_WARN   -- may not be wrong, but undesirable stuff happens.  *
 * - MJOLNIR_LOG_ERROR  -- something wrong.                                  *
 *
 * - MJOLNIR_SCOPE      -- helps logging by indentation.                     *
 * - MJOLNIR_GET_LOGGER -- get a logger for the current scope with name.     *
 * - MJOLNIR_GET_DEFAULT_LOGGER -- get a default logger.                     */

namespace logger_detail
{
// since there is no standard way to format containers, the standard library
// does not provide output operators for containers. To make logging easier,
// here output operators are defined for some of the containers.
// To avoid fixing the output format, these are defined in the special namespace
// `logger_detail`. loggers first imports this namespace and output using these
// operators.

template<typename charT, typename traits, typename T, std::size_t N>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const std::array<T, N>& ar)
{
    os << '[';
    for(const auto& v : ar){os << v << ", ";}
    os << ']';
    return os;
}

template<typename charT, typename traits, typename T, typename Alloc>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os,
           const std::vector<T, Alloc>& vc)
{
    os << '[';
    for(const auto item : vc){os << item << ", ";}
    os << ']';
    return os;
}

template<typename charT, typename traits, typename Key, typename Value,
         typename Comp, typename Alloc>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os,
           const std::map<Key, Value, Comp, Alloc>& mp)
{
    os << '{';
    for(const auto kv : mp){os << kv.first << '=' << kv.second << ", ";}
    os << '}';
    return os;
}

template<typename charT, typename traits, typename Key, typename Value,
         typename Hash, typename Pred, typename Alloc>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os,
           const std::unordered_map<Key, Value, Hash, Pred, Alloc>& mp)
{
    os << '{';
    for(const auto kv : mp){os << kv.first << '=' << kv.second << ", ";}
    os << '}';
    return os;
}

template<typename charT, typename traits, typename T1, typename T2>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const std::pair<T1, T2>& pr)
{
    os << '[' << pr.first << ", " << pr.second << ']';
    return os;
}

template<std::size_t I, std::size_t N>
struct tuple_output_helper
{
    template<typename charT, typename traits, typename ... Ts>
    static std::basic_ostream<charT, traits>&
    invoke(std::basic_ostream<charT, traits>& os, const std::tuple<Ts...>& t)
    {
        static_assert(N == sizeof...(Ts), "");
        os << std::get<I>(t) << ", ";
        return tuple_output_helper<I+1, N>::invoke(os, t);
    }
};
template<std::size_t N>
struct tuple_output_helper<N, N>
{
    template<typename charT, typename traits, typename ... Ts>
    static std::basic_ostream<charT, traits>&
    invoke(std::basic_ostream<charT, traits>& os, const std::tuple<Ts...>& t)
    {
        static_assert(N == sizeof...(Ts), "");
        return os;
    }
};
template<typename charT, typename traits, typename ... Ts>
std::basic_ostream<charT, traits>&
operator<<(std::basic_ostream<charT, traits>& os, const std::tuple<Ts...>& t)
{
    os << '[';
    tuple_output_helper<0, sizeof...(Ts)>::invoke(os, t);
    os << ']';
    return os;
}

} // logger_detail

template<typename charT, typename traits = std::char_traits<charT>>
class basic_logger
{
  public:

    using string_type  = std::basic_string <charT, traits>;
    using fstream_type = std::basic_fstream<charT, traits>;
    using ostream_type = std::basic_ostream<charT, traits>;

    static constexpr std::size_t indent_size = 2;

    enum class Level : std::uint8_t
    {
        None   = 0, // internal use only
        Debug  = 1,
        Info   = 2,
        Notice = 3,
        Warn   = 4,
        Error  = 5
    };

  public:

    explicit basic_logger(const std::string& fname)
        : indentation_needed_(false), indent_(0), fname_(fname)
    {
        // clear the contents and check the file can be opened
        fstream_type ofs(this->fname_, std::ios_base::out);
        if(!ofs.good())
        {
            std::cerr << "logger: file open error: " << fname_ << std::endl;
            std::exit(1);
        }
    }
    ~basic_logger() = default;

    void indent()   noexcept {indent_ += 1; return;}
    void unindent() noexcept {indent_ -= 1; return;}

    template<typename ... Ts>
    void log(Level level, Ts&& ... args)
    {
        fstream_type ofs(this->fname_, std::ios_base::out | std::ios_base::app);
        if(!ofs.good())
        {
            throw_exception<std::runtime_error>(
                "Logger: file open error: ", this->fname_);
        }
        if(this->indentation_needed_)
        {
            ofs << string_type(indent_size * indent_, ' ');
        }

        if(level == Level::Notice || level == Level::Warn || level == Level::Error)
        {
            // message is also printed to stderr
            output_message(std::cerr, "-- ", this->stringize(level),
                           std::forward<Ts>(args)...);
            std::cerr << std::endl;
        }

        output_message(ofs, this->stringize(level), std::forward<Ts>(args)...);
        ofs << std::endl;

        this->indentation_needed_ = true;
        return;
    }

    template<typename ... Ts>
    void log_no_lf(Level level, Ts&& ... args)
    {
        fstream_type ofs(this->fname_, std::ios_base::out | std::ios_base::app);
        if(!ofs.good())
        {
            throw_exception<std::runtime_error>(
                "Logger: file open error: ", this->fname_);
        }
        if(this->indentation_needed_)
        {
            ofs << string_type(indent_size * indent_, ' ');
        }

        if(level == Level::Notice || level == Level::Warn || level == Level::Error)
        {
            // message is also printed to stderr
            output_message(std::cerr, "-- ", this->stringize(level),
                           std::forward<Ts>(args)...);
            std::cerr << std::endl;
        }

        output_message(ofs, this->stringize(level), std::forward<Ts>(args)...);
        ofs << std::flush;

        this->indentation_needed_ = false;
        return;
    }

  private:

    std::string stringize(const Level lv)
    {
        switch(lv)
        {
            case Level::None  : return "";
            case Level::Debug : return "";
            case Level::Info  : return "";
            case Level::Notice: return "";
            case Level::Warn  : return "WARNING: ";
            case Level::Error : return "ERROR: ";
            default: return "UNKNOWN LOG LEVEL: ";
        }
    }

    template<typename T, typename ...T_args>
    static void output_message(ostream_type& os, T&& arg1, T_args&& ...args)
    {
        using namespace logger_detail; // to output containers
        os << arg1;
        return output_message(os, std::forward<T_args>(args)...);
    }
    static void output_message(ostream_type& os) {return;}

  private:

    bool indentation_needed_;
    std::size_t indent_;
    string_type  fname_;
};

template<typename charT, typename traitsT = std::char_traits<charT>>
class basic_logger_manager
{
  public:
    typedef basic_logger<charT, traitsT>         logger_type;
    typedef std::unique_ptr<logger_type>         resource_type;
    typedef std::map<std::string, resource_type> container_type;

  public:

    static void set_default_logger(const std::string& fname)
    {
        if(default_ == fname)
        {
            std::cerr << "WARNING for Developers: Default Logger("
                      << fname << ") is set twice." << std::endl;
            return;
        }

        default_ = fname;
        if(loggers_.count(fname) != 0)
        {
            std::cerr << "WARNING for Developers: Logger(" << fname << ") is "
                      << "already set. from now, it becomes the default logger."
                      << std::endl;
            return;
        }
        loggers_.emplace(fname, make_unique<logger_type>(fname));
        return;
    }

    static logger_type& get_default_logger()
    {
        if(default_.empty())
        {
            throw_exception<std::out_of_range>("mjolnir::basic_logger_manager: "
                "default logger is not set yet.");
        }
        if(loggers_.count(default_) == 0)
        {
            throw_exception<std::out_of_range>("mjolnir::basic_logger_manager: "
                "default logger (", default_, ") does not exist");
        }
        return *(loggers_.at(default_));
    }

    static logger_type& get_logger(const std::string& name)
    {
        if(loggers_.count(name) == 0)
        {
            loggers_.emplace(name, make_unique<logger_type>(name));
        }
        return *(loggers_.at(name));
    }

  private:

    static std::string    default_;
    static container_type loggers_;
};

template<typename charT, typename traitsT>
typename basic_logger_manager<charT, traitsT>::container_type
basic_logger_manager<charT, traitsT>::loggers_;

template<typename charT, typename traitsT>
std::string basic_logger_manager<charT, traitsT>::default_;

template<typename charT, typename traitsT = std::char_traits<charT>>
class basic_scope
{
  public:
    typedef charT   char_type;
    typedef traitsT traits_type;
    typedef basic_logger<char_type, traits_type> logger_type;

  public:

    basic_scope(logger_type& trc, const std::string& name)
      : start_(std::chrono::system_clock::now()), logger_(trc), name_(name)
    {
        logger_.log(logger_type::Level::None, true, this->name_, " {\n");
        logger_.indent();
    }
    ~basic_scope()
    {
        logger_.unindent();
        logger_.log(logger_type::Level::None, true, "} ", this->format_duration(
                    std::chrono::system_clock::now() - this->start_), '\n');
    }

    std::string const& name() const noexcept {return name_;}

  private:

    std::string format_duration(const std::chrono::system_clock::duration& dur)
    {
        const auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(dur);
        if(ns.count() < 1000)
        {
            return std::to_string(ns.count()) + " [ns]";
        }
        else if(ns.count() < 1000000) // 10^6, less than milliseconds
        {
            return std::to_string(ns.count() * 1e-3) + " [us]";
        }
        else if(ns.count() < 1000000000) // 10^9, less than seconds
        {
            return std::to_string(ns.count() * 1e-6) + " [ms]";
        }
        else
        {
            return std::to_string(ns.count() * 1e-9) + " [sec]";
        }
    }

  private:
    std::chrono::system_clock::time_point start_;
    logger_type& logger_;
    std::string name_;
};

using Logger        = basic_logger<char>;
using LoggerManager = basic_logger_manager<char>;
using Scope         = basic_scope<char>;

#define MJOLNIR_SET_DEFAULT_LOGGER(name)  LoggerManager::set_default_logger(name)
#define MJOLNIR_GET_DEFAULT_LOGGER()      auto& l_o_g_g_e_r_ = LoggerManager::get_default_logger()
#define MJOLNIR_GET_LOGGER(name)          auto& l_o_g_g_e_r_ = LoggerManager::get_logger(name)
#define MJOLNIR_LOG_INFO(args...)         l_o_g_g_e_r_.log(Logger::Level::Info,   args)
#define MJOLNIR_LOG_NOTICE(args...)       l_o_g_g_e_r_.log(Logger::Level::Notice, args)
#define MJOLNIR_LOG_WARN(args...)         l_o_g_g_e_r_.log(Logger::Level::Warn,   args)
#define MJOLNIR_LOG_ERROR(args...)        l_o_g_g_e_r_.log(Logger::Level::Error,  args)
#define MJOLNIR_LOG_INFO_NO_LF(args...)   l_o_g_g_e_r_.log_no_lf(Logger::Level::Info,   args)
#define MJOLNIR_LOG_NOTICE_NO_LF(args...) l_o_g_g_e_r_.log_no_lf(Logger::Level::Notice, args)
#define MJOLNIR_LOG_WARN_NO_LF(args...)   l_o_g_g_e_r_.log_no_lf(Logger::Level::Warn,   args)
#define MJOLNIR_LOG_ERROR_NO_LF(args...)  l_o_g_g_e_r_.log_no_lf(Logger::Level::Error,  args)
#define MJOLNIR_SCOPE(name, id)           Scope s_c_o_p_e_##id (l_o_g_g_e_r_, #name)

#ifdef MJOLNIR_DEBUG
#  define MJOLNIR_GET_DEFAULT_LOGGER_DEBUG() auto& l_o_g_g_e_r_ = LoggerManager::get_default_logger()
#  define MJOLNIR_LOG_DEBUG(args...)         l_o_g_g_e_r_.log(Logger::Level::Debug, args)
#  define MJOLNIR_LOG_DEBUG_NO_LF(args...)   l_o_g_g_e_r_.log_no_lf(Logger::Level::Debug, args)
#  define MJOLNIR_SCOPE_DEBUG(name, id)      Scope s_c_o_p_e_##id (l_o_g_g_e_r_, #name)
#else
#  define MJOLNIR_GET_DEFAULT_LOGGER_DEBUG() /**/
#  define MJOLNIR_LOG_DEBUG(args...)         /**/
#  define MJOLNIR_LOG_DEBUG_NO_LF(args...)   /**/
#  define MJOLNIR_SCOPE_DEBUG(name, id)      /**/
#endif // MJOLNIR_DEBUG

} // mjolnir
#endif /* MJOLNIR_UTIL_LOGGER */
