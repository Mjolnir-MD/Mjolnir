#ifndef MJOLNIR_UTIL_LOGGER
#define MJOLNIR_UTIL_LOGGER
#include "make_unique.hpp"
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <chrono>
#include <ctime>

namespace mjolnir
{

template<typename charT, typename traits = std::char_traits<charT>,
         typename alloc = std::allocator<charT>>
class basic_logger
{
  public:
    typedef std::basic_string<charT, traits, alloc>        string_type;
    typedef std::basic_ostream<charT, traits>              ostream_type;
    typedef std::basic_fstream<charT, traits>              fstream_type;
    typedef std::basic_ostringstream<charT, traits, alloc> ostringstream_type;
    constexpr static std::size_t name_length = 25;

    enum class Level
    {
        Debug,
        Time,
        Info,
        Warn,
        Error,
    };

  public:

    explicit basic_logger(const string_type& name)
        : name_(name.size() <= name_length ? name : name.substr(0, name_length))
    {}
    ~basic_logger() = default;

    template<typename ... T_args>
    void log(Level level, T_args&& ... args) const
    {
        if(output_.count(level) == 1)
        {
            fstream_type ofs(output_.at(level),
                             fstream_type::out | fstream_type::app);
            if(!ofs.good())
                throw std::runtime_error(
                        "file open erorr: " + output_.at(level));

            log_(ofs, level, std::forward<T_args>(args)...);
            ofs.close();
            return;
        }
        else
        {
            log_(std::cerr, level, std::forward<T_args>(args)...);
            return;
        }
    }

    void set_output(Level lev, std::string&& fname)
    {
        output_[lev] = std::forward<string_type>(fname);
        return;
    }

  private:

    template<typename ... T_args>
    void log_(ostream_type& os, const Level& level, T_args&& ...args) const
    {
        std::time_t t = std::chrono::system_clock::to_time_t(
                std::chrono::system_clock::now());
        std::tm tm = *std::localtime(&t);

        os << to_string(level) << '|' << to_str(tm) << "| ";
        os << std::setw(name_length) << std::left << name_;
        string_type mes = gen_message(std::forward<T_args>(args)...);
        os << mes << std::endl;
        return;
    }

    template<typename T, typename ...T_args>
    string_type gen_message(T&& arg1, T_args&& ...args) const
    {
        ostringstream_type oss;
        oss << " " << arg1;
        return oss.str() + gen_message(std::forward<T_args>(args)...);
    }

    template<typename T>
    string_type gen_message(T&& arg1) const
    {
        ostringstream_type oss;
        oss << " " << arg1;
        return oss.str();
    }

    string_type to_string(Level l) const
    {
        switch(l)
        {
            case Level::Debug: return "[Debug  ]";
            case Level::Time:  return "[Time   ]";
            case Level::Info:  return "[Info   ]";
            case Level::Warn:  return "[Warning]";
            case Level::Error: return "[Error  ]";
            default:           return "[Unknown]";
        }
    }

    string_type to_str(const std::tm& t) const
    {
        ostringstream_type oss;
        oss << t.tm_year+1900 << '/' << t.tm_mon+1 << '/' << t.tm_mday << ' ';
        oss << t.tm_hour << ':' << t.tm_min << ':' << t.tm_sec;
        return oss.str();
    }

  private:

    const string_type name_;
    std::map<Level, std::string> output_;
};

template<typename charT, typename char_traits = std::char_traits<charT>,
         typename alloc = std::allocator<charT>>
class LoggerManager
{
  public:
    typedef std::basic_string<charT, char_traits, alloc> string_type;
    typedef basic_logger<charT, char_traits, alloc>      logger_type;
    typedef std::unique_ptr<logger_type>                 resource_type;
    typedef std::map<string_type, resource_type>         container_type;
    constexpr static const char* debug_file = "mjolnir_debug.log";
    constexpr static const char* info_file  = "mjolnir_info.log";

  public:
    static logger_type& get_logger(const string_type& name)
    {
        if(dirty)
        {
            std::ofstream dbg(debug_file);
            std::ofstream info(info_file);
            dbg.close();
            info.close();
            dirty = false;
        }
        if(loggers_.count(name) == 0)
        {
            resource_type newlogger = make_unique<logger_type>(name);
            newlogger->set_output(logger_type::Level::Debug, debug_file);
            newlogger->set_output(logger_type::Level::Info,  info_file);
            loggers_.emplace(name, std::move(newlogger));
        }
        return *(loggers_.at(name));
    }

  private:
    static bool dirty;
    static container_type loggers_;
};

template<typename charT, typename traits, typename alloc>
typename LoggerManager<charT, traits, alloc>::container_type
LoggerManager<charT, traits, alloc>::loggers_;

template<typename charT, typename traits, typename alloc>
bool LoggerManager<charT, traits, alloc>::dirty = true;

// duration_as_string {{{
template<typename durationT>
struct duration_as_string;

template<>
struct duration_as_string<std::chrono::nanoseconds>
{
    static std::string invoke(){return "[nsec]";}
};

template<>
struct duration_as_string<std::chrono::microseconds>
{
    static std::string invoke(){return "[usec]";}
};

template<>
struct duration_as_string<std::chrono::milliseconds>
{
    static std::string invoke(){return "[msec]";}
};

template<>
struct duration_as_string<std::chrono::seconds>
{
    static std::string invoke(){return "[sec ]";}
};

template<>
struct duration_as_string<std::chrono::minutes>
{
    static std::string invoke(){return "[min ]";}
};

template<>
struct duration_as_string<std::chrono::hours>
{
    static std::string invoke(){return "[hour]";}
};
// }}}

template<typename durationT,
         typename charT = char, typename traits = std::char_traits<charT>,
         typename alloc = std::allocator<charT>>
struct stopwatch
{
    typedef duration_as_string<durationT> as_string;
    typedef basic_logger<charT, traits, alloc> logger_type;

    stopwatch(const std::string& name, logger_type& logger)
        : start(std::chrono::system_clock::now()), name_(name), logger_(logger)
    {}
    ~stopwatch()
    {
        const auto stop = std::chrono::system_clock::now();
        const auto dur = std::chrono::duration_cast<durationT>(stop - start);
        logger_.log(logger_type::Time, name_, dur.count(), as_string::invoke());
    }

  private:

    const std::chrono::system_clock::time_point start;
    const std::string name_;
    logger_type& logger_;
};

typedef basic_logger<char>     Logger;
typedef basic_logger<wchar_t>  wLogger;
typedef basic_logger<char16_t> u16Logger;
typedef basic_logger<char32_t> u32Logger;

typedef basic_logger<char>::Level     LogLevel;
typedef basic_logger<wchar_t>::Level  wLogLevel;
typedef basic_logger<char16_t>::Level u16LogLevel;
typedef basic_logger<char32_t>::Level u32LogLevel;

#ifdef MJOLNIR_DEBUG
#define MJOLNIR_SET_LOGGER(name)   auto& logger_ = LoggerManager<char>::get_logger(name);
#define MJOLNIR_LOG(args...)       logger_.log(args);
#define MJOLNIR_LOG_DEBUG(args...) logger_.log(LogLevel::Debug, args);
#define MJOLNIR_LOG_INFO(args...)  logger_.log(LogLevel::Info,  args);
#define MJOLNIR_LOG_WARN(args...)  logger_.log(LogLevel::Warn,  args);
#define MJOLNIR_LOG_ERROR(args...) logger_.log(LogLevel::Error, args);
#elif defined(MJOLNIR_NO_DEBUG)
#define MJOLNIR_SET_LOGGER(name)   /*name*/
#define MJOLNIR_LOG(args...)       /*args*/
#define MJOLNIR_LOG_DEBUG(args...) /*args*/
#define MJOLNIR_LOG_INFO(args...)  /*args*/
#define MJOLNIR_LOG_WARN(args...)  /*args*/
#define MJOLNIR_LOG_ERROR(args...) /*args*/
#else  // normal one
#define MJOLNIR_SET_LOGGER(name)   auto& logger_ = LoggerManager<char>::get_logger(name);
#define MJOLNIR_LOG(args...)       /*args*/
#define MJOLNIR_LOG_DEBUG(args...) /*args*/
#define MJOLNIR_LOG_INFO(args...)  logger_.log(LogLevel::Info,  args);
#define MJOLNIR_LOG_WARN(args...)  logger_.log(LogLevel::Warn,  args);
#define MJOLNIR_LOG_ERROR(args...) logger_.log(LogLevel::Error, args);
#endif // MJOLNIR_DEBUG

} // mjolnir
#endif /* MJOLNIR_UTIL_LOGGER */
