#ifndef MJOLNIR_UTIL_LOGGER
#define MJOLNIR_UTIL_LOGGER
#include <mjolnir/util/make_unique.hpp>
#include <mjolnir/util/throw_exception.hpp>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <chrono>

namespace mjolnir
{

template<typename charT, typename traits = std::char_traits<charT>>
class basic_logger
{
  public:

    using string_type  = std::basic_string <charT, traits>;
    using fstream_type = std::basic_fstream<charT, traits>;
    using ostream_type = std::basic_ostream<charT, traits>;

    static constexpr std::size_t indent_size = 2;
    enum class Level
    {
        Debug,
        Warn,
        Error
    };

  public:

    explicit basic_logger(const std::string& fname): indent_(0), fname_(fname)
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
    void log(Level level, Ts&& ... args) const
    {
        fstream_type ofs(this->fname_, std::ios_base::out | std::ios_base::app);
        if(!ofs.good())
        {
            throw_exception<std::runtime_error>("Logger: file open error: ",
                    this->fname_);
        }
        ofs << string_type(indent_size * indent_, ' ');
        if     (level == Level::Warn)  {ofs << "[WARNING]";}
        else if(level == Level::Error) {ofs << "[ERROR]";}
        output_message(ofs, std::forward<Ts>(args)...);
        ofs << std::endl;
        return;
    }

    template<typename ... Ts>
    void write(Ts&& ... args) const
    {
        fstream_type ofs(this->fname_, std::ios_base::out | std::ios_base::app);
        if(!ofs.good())
        {
            throw_exception<std::runtime_error>("Logger: file open error: ",
                    this->fname_);
        }
        ofs << string_type(indent_size * indent_, ' ');
        output_message(ofs, std::forward<Ts>(args)...);
        ofs << std::endl;
        return;
    }

  private:

    template<typename T, typename ...T_args>
    static void output_message(ostream_type& os, T&& arg1, T_args&& ...args)
    {
        os << ' ' << arg1;
        return output_message(os, std::forward<T_args>(args)...);
    }
    static void output_message(ostream_type& os) {return;}

  private:

    std::size_t indent_;
    string_type  fname_;
};

template<typename charT, typename traitsT = std::char_traits<charT>>
class basic_logger_manager
{
  public:
    typedef basic_logger<charT, traitsT>     logger_type;
    typedef std::unique_ptr<logger_type>         resource_type;
    typedef std::map<std::string, resource_type> container_type;

  public:

    static logger_type& get_logger(const std::string name)
    {
        if(loggers_.count(name) == 0)
        {
            loggers_.emplace(name, make_unique<logger_type>(name));
        }
        return *(loggers_.at(name));
    }

  private:

    static container_type loggers_;
};

template<typename charT, typename traitsT>
typename basic_logger_manager<charT, traitsT>::container_type
basic_logger_manager<charT, traitsT>::loggers_;

template<typename charT, typename traitsT = std::char_traits<charT>>
class basic_scope
{
  public:
    typedef charT   char_type;
    typedef traitsT traits_type;
    typedef basic_logger<char_type, traits_type> logger_type;

  public:

    basic_scope(logger_type& trc, const std::string& name)
      : logger_(trc), name_(name)
    {
        logger_.write(name_, " {");
        logger_.indent();
    }
    ~basic_scope()
    {
        logger_.unindent();
        logger_.write('}');
    }

    std::string const& name() const noexcept {return name_;}

  private:
    logger_type& logger_;
    std::string name_;
};

using Logger        = basic_logger<char>;
using LoggerManager = basic_logger_manager<char>;
using Scope         = basic_scope<char>;

#ifdef MJOLNIR_DEBUG
#  define MJOLNIR_SET_DEFAULT_LOGGER()\
      auto& l_o_g_g_e_r_ = LoggerManager::get_logger("mjolnir.log")
#  define MJOLNIR_SET_LOGGER(name)\
      auto& l_o_g_g_e_r_ = LoggerManager::get_logger(name)
#  define MJOLNIR_SCOPE(name, id)    Scope s_c_o_p_e_##id (l_o_g_g_e_r_, #name)
#  define MJOLNIR_LOG_DEBUG(args...) l_o_g_g_e_r_.log(Logger::Level::Debug, args)
#  define MJOLNIR_LOG_WARN(args...)  l_o_g_g_e_r_.log(Logger::Level::Warn,  args)
#  define MJOLNIR_LOG_ERROR(args...) l_o_g_g_e_r_.log(Logger::Level::Error, args)
#else
#  define MJOLNIR_SET_DEFAULT_LOGGER() /**/
#  define MJOLNIR_SET_LOGGER(name)     /**/
#  define MJOLNIR_SCOPE(name, id)      /**/
#  define MJOLNIR_LOG_DEBUG(args...)   /**/
#  define MJOLNIR_LOG_WARN(args...)    /**/
#  define MJOLNIR_LOG_ERROR(args...)   /**/
#endif // MJOLNIR_DEBUG

} // mjolnir
#endif /* MJOLNIR_UTIL_LOGGER */
