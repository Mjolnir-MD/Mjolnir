#ifndef MJOLNIR_UTIL_LOGGER
#define MJOLNIR_UTIL_LOGGER
#include "make_unique.hpp"
#include <string>
#include <fstream>
#include <iostream>
#include <map>

namespace mjolnir
{

template<typename charT, typename char_traits = std::char_traits<charT>>
class basic_logger
{
  public:
    typedef std::basic_string<charT, char_traits> string_type;
    typedef std::basic_ostream<charT, char_traits> ostream_type;
    typedef std::basic_fstream<charT, char_traits> fstream_type;

    enum class Level
    {
        Debug,
        Info,
        Warn,
        Error,
    };

  public:

    explicit basic_logger(const string_type& name): name_(name){}
    basic_logger(basic_logger&& l) = default;
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

            log_(ofs, to_string(level), name_, std::forward<T_args>(args)...);
            ofs.close();
            return;
        }
        else
        {// TODO wchar, u16char, u32char
            return log_(std::cerr, to_string(level), name_,
                        std::forward<T_args>(args)...);
        }
    }

    void set_output(Level lev, string_type&& fname)
    {
        output_[lev] = std::forward<string_type>(fname);
        return;
    }

  private:

    template<typename T, typename ... T_args>
    void log_(ostream_type& os, T&& arg1, T_args&& ...args) const
    {
        os << std::forward<T>(arg1) << " ";
        log_(os, args...);
        return;
    }

    template<typename T>
    void log_(ostream_type& os, T&& arg1) const
    {
        os << std::forward<T>(arg1) << std::endl;
        return;
    }

    string_type to_string(Level l) const
    {
        switch(l)
        {
            case Level::Debug: return "Log-Debug:";
            case Level::Info:  return "Log-Info:";
            case Level::Warn:  return "Log-Warn:";
            case Level::Error: return "Log-Error:";
            default:           return "Log-Unknown";
        }
    }

  private:

    string_type name_;
    std::map<Level, string_type> output_;

};

template<typename charT, typename char_traits = std::char_traits<charT>>
class LoggerManager
{
  public:
    typedef std::basic_string<charT, char_traits> string_type;

  public:
    static
    basic_logger<charT, char_traits>& get_logger(const string_type& name)
    {
        if(loggers_.count(name) == 0)
            loggers_.emplace(
                    name, make_unique<basic_logger<charT, char_traits>>(name));
        return *(loggers_[name]);
    }

  private:
    static
    std::map<string_type,
             std::unique_ptr<basic_logger<charT, char_traits>>> loggers_;
};

template<typename charT, typename traits>
std::map<std::basic_string<charT, traits>,
         std::unique_ptr<basic_logger<charT, traits>>>
LoggerManager<charT, traits>::loggers_;

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
