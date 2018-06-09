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

/* Here, these macros are defined.                                   *
 *  1. MJOLNIR_LOG_DEBUG                                             *
 *    - just for debug purpose. the log file may become too large.   *
 *  2. MJOLNIR_LOG_INFO                                              *
 *    - useful information to check the behavior. by default, on.    *
 *  3. MJOLNIR_LOG_WARN                                              *
 *    - be careful with it. It can run, but it might not make sense. *
 *  4. MJOLNIR_LOG_ERROR                                             *
 *    - something wrong happen.                                      *
 *  5. MJOLNIR_SCOPE                                                 *
 *    - helps output format of log.                                  *
 *  6. MJOLNIR_SCOPE_DEBUG                                           *
 *    - helps output format of log. by default, off.                 *
 *  7. MJOLNIR_GET_DEFAULT_LOGGER                                    *
 *    - enable logging in the scope.                                 *
 *  8. MJOLNIR_GET_LOGGER                                            *
 *    - make new log file and output to that file.                   *
 *  9. MJOLNIR_GET_DEFAULT_LOGGER_DEBUG                              *
 *    - enabled only when MJOLNIR_DEBUG flag is set.                 *
 * 10. MJOLNIR_GET_LOGGER_DEBUG                                      *
 *    - enabled only when MJOLNIR_DEBUG flag is set.                 */
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
        Info,
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
        if(indent_ != 0)
        {
            ofs << string_type(indent_size * indent_, ' ');
        }
        if     (level == Level::Warn)  {ofs << "[WARNING] ";}
        else if(level == Level::Error) {ofs << "[ERROR] ";}
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
        if(indent_ != 0)
        {
            ofs << string_type(indent_size * indent_, ' ');
        }
        output_message(ofs, std::forward<Ts>(args)...);
        ofs << std::endl;
        return;
    }

  private:

    template<typename T, typename ...T_args>
    static void output_message(ostream_type& os, T&& arg1, T_args&& ...args)
    {
        os << arg1;
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
    typedef basic_logger<charT, traitsT>         logger_type;
    typedef std::unique_ptr<logger_type>         resource_type;
    typedef std::map<std::string, resource_type> container_type;

  public:

    static void set_default_logger(const std::string& fname)
    {
        default_ = fname;
        if(loggers_.count(fname) == 0)
        {
            loggers_.emplace(fname, make_unique<logger_type>(fname));
            return;
        }
        else
        {
            throw_exception<std::logic_error>("mjolnir::basic_logger_manager: "
                "default logger `", fname, "` is already set");
        }
    }

    static logger_type& get_default_logger()
    {
        if(loggers_.count(default_) == 0)
        {
            throw_exception<std::out_of_range>("mjolnir::basic_logger_manager: "
                "default logger `", default_, "` does not exist");
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
#  define MJOLNIR_SET_DEFAULT_LOGGER(name)\
    LoggerManager::set_default_logger(name)
#  define MJOLNIR_GET_DEFAULT_LOGGER()\
    auto& l_o_g_g_e_r_ = LoggerManager::get_default_logger()
#  define MJOLNIR_GET_DEFAULT_LOGGER_DEBUG()\
    auto& l_o_g_g_e_r_ = LoggerManager::get_default_logger()
#  define MJOLNIR_GET_LOGGER(name)\
    auto& l_o_g_g_e_r_ = LoggerManager::get_logger(name)
#  define MJOLNIR_SCOPE(name, id)       Scope s_c_o_p_e_##id (l_o_g_g_e_r_, #name)
#  define MJOLNIR_SCOPE_DEBUG(name, id) Scope s_c_o_p_e_##id (l_o_g_g_e_r_, #name)
#  define MJOLNIR_LOG_DEBUG(args...) l_o_g_g_e_r_.log(Logger::Level::Debug, args)
#  define MJOLNIR_LOG_INFO(args...)  l_o_g_g_e_r_.log(Logger::Level::Info,  args)
#  define MJOLNIR_LOG_WARN(args...)  l_o_g_g_e_r_.log(Logger::Level::Warn,  args)
#  define MJOLNIR_LOG_ERROR(args...) l_o_g_g_e_r_.log(Logger::Level::Error, args)
#elif defined(MJOLNIR_NO_LOG)
#  define MJOLNIR_SET_DEFAULT_LOGGER(name)   /**/
#  define MJOLNIR_GET_DEFAULT_LOGGER()       /**/
#  define MJOLNIR_GET_DEFAULT_LOGGER_DEBUG() /**/
#  define MJOLNIR_GET_LOGGER(name)           /**/
#  define MJOLNIR_SCOPE(name, id)            /**/
#  define MJOLNIR_SCOPE_DEBUG(name, id)      /**/
#  define MJOLNIR_LOG_DEBUG(args...)         /**/
#  define MJOLNIR_LOG_INFO(args...)          /**/
#  define MJOLNIR_LOG_WARN(args...)          /**/
#  define MJOLNIR_LOG_ERROR(args...)         /**/
#else
#  define MJOLNIR_SET_DEFAULT_LOGGER(name)\
    LoggerManager::set_default_logger(name)
#  define MJOLNIR_GET_DEFAULT_LOGGER()\
    auto& l_o_g_g_e_r_ = LoggerManager::get_default_logger()
#  define MJOLNIR_GET_LOGGER(name)\
    auto& l_o_g_g_e_r_ = LoggerManager::get_logger(name)
#  define MJOLNIR_SCOPE(name, id) Scope s_c_o_p_e_##id (l_o_g_g_e_r_, #name)
#  define MJOLNIR_GET_DEFAULT_LOGGER_DEBUG() /**/
#  define MJOLNIR_SCOPE_DEBUG(name, id)      /**/
#  define MJOLNIR_LOG_DEBUG(args...)         /**/
#  define MJOLNIR_LOG_INFO(args...)  l_o_g_g_e_r_.log(Logger::Level::Info,  args)
#  define MJOLNIR_LOG_WARN(args...)  l_o_g_g_e_r_.log(Logger::Level::Warn,  args)
#  define MJOLNIR_LOG_ERROR(args...) l_o_g_g_e_r_.log(Logger::Level::Error, args)
#endif // MJOLNIR_DEBUG

} // mjolnir
#endif /* MJOLNIR_UTIL_LOGGER */
