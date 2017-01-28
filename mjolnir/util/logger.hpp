#ifndef MJOLNIR_UTIL_LOGGER
#define MJOLNIR_UTIL_LOGGER
#include <string>
#include <fstream>
#include <iostream>
#include <map>
#include "make_unique.hpp"

namespace mjolnir
{

template<typename charT, typename char_traits = std::char_traits<charT>>
class Logger
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

    explicit Logger(const string_type& name): name_(name){}
    Logger(Logger&& l) = default;
    ~Logger() = default;

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
        {
            return log_(std::cerr, to_string(level), name_,
                        std::forward<T_args>(args)...);
        }
    }

    void set_output(Level lev, string_type&& fname)
    {
        output_[lev] = std::forward<string_type>(fname);
        return;
    }

    static
    Logger& get_logger(const string_type& name)
    {
        if(loggers_.count(name) == 0)
            loggers_.emplace(name, make_unique<Logger>(name));
        return *(loggers_[name]);
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
            case Level::Debug: return "Log-Debug: ";
            case Level::Info:  return "Log-Info:  ";
            case Level::Warn:  return "Log-Warn:  ";
            case Level::Error: return "Log-Error: ";
            default:           return "Log-Unknown ";
        }
    }

  private:

    string_type name_;
    std::map<Level, string_type> output_;

    static
    std::map<string_type, std::unique_ptr<Logger>> loggers_;
};

template<typename charT, typename traits>
std::map<std::basic_string<charT, traits>,
         std::unique_ptr<Logger<charT, traits>>>
Logger<charT, traits>::loggers_;

#ifdef MJOLNIR_DEBUG
#define MJOLNIR_LOG(args...) logger_.log(args);
#else  // NO-DEBUG
#define MJOLNIR_LOG(args...) /* args */
#endif // MJOLNIR_DEBUG

} // mjolnir
#endif /* MJOLNIR_UTIL_LOGGER */
