#ifndef MJOLNIR_UTIL_STATIC_STRING_HPP
#define MJOLNIR_UTIL_STATIC_STRING_HPP
#include <mjolnir/util/throw_exception.hpp>
#include <string>
#include <stdexcept>
#include <istream>
#include <ostream>
#include <cstdint>

namespace mjolnir
{

// static_string is a statically sized std::string, inspired by static_vector
// provided by Boost.Container.
// this is the minimal implementation for a specific perpose, so it doesn't have
// all the functionalities that exists in `std::string`. For example, it assumes
// that the character encoding is ascii. Basically, it can contain UTF-8 string
// as byte-array. But in that case, `size()` might not fit to the actual number
// of letters.
template<std::size_t N>
class static_string
{
    // because it contains '\0' delimiter at the end of buffer, the actual
    // free buffer size is N-1. static_string is to manage the size of object,
    // the clarity of the size is considered more important than the usability.
    static constexpr std::size_t capacity_ = N - 1;

    // to make code shorter.
    void check_size() const
    {
        if(this->size_ > this->capacity_)
        {
            throw_exception<std::length_error>("static_string (with capacity ",
                this->capacity_, ") cannot be constructed with string having ",
                this->size_, " length.");
        }
        return ;
    }

  public:

    using char_type              = char;
    using traits_type            = std::char_traits<char>;
    using container_type         = std::array<char, N>;

    using value_type             = typename container_type::value_type;
    using size_type              = typename container_type::size_type;
    using difference_type        = typename container_type::difference_type;
    using pointer                = typename container_type::pointer;
    using const_pointer          = typename container_type::const_pointer;
    using reference              = typename container_type::reference;
    using const_reference        = typename container_type::const_reference;
    using iterator               = typename container_type::iterator;
    using const_iterator         = typename container_type::const_iterator ;
    using reverse_iterator       = typename container_type::reverse_iterator;
    using const_reverse_iterator = typename container_type::const_reverse_iterator;

  public:

    static_string() : size_(0) {buffer_.fill('\0');}
    ~static_string() = default;
    static_string(static_string const&) = default;
    static_string(static_string &&)     = default;
    static_string& operator=(static_string const&) = default;
    static_string& operator=(static_string &&)     = default;

    explicit static_string(const char* literal): size_(traits_type::length(literal))
    {
        this->check_size();
        traits_type::copy(this->data(), literal, this->size_+1);
    }
    explicit static_string(const std::string& str): size_(str.size())
    {
        this->check_size();
        traits_type::copy(this->data(), str.c_str(), this->size_+1);
    }
    template<std::size_t M>
    explicit static_string(const static_string<M>& str): size_(str.size())
    {
        if(N < M) {this->check_size();}
        traits_type::copy(this->data(), str.c_str(), this->size_+1);
    }

    static_string& operator=(const char* literal)
    {
        this->size_ = traits_type::length(literal);
        this->check_size();
        traits_type::copy(this->data(), literal, this->size_+1);
        return *this;
    }
    static_string& operator=(const std::string& str)
    {
        this->size_ = str.size();
        this->check_size();
        traits_type::copy(this->data(), str.c_str(), this->size_+1);
        return *this;
    }
    template<std::size_t M>
    static_string& operator=(static_string<M> const& rhs)
    {
        this->size_ = rhs.size();
        if(N < M) {this->check_size();}
        traits_type::copy(this->data(), rhs.c_str(), this->size_+1);
        return *this;
    }

    static_string& operator+=(const char c)
    {
        this->size_++;
        this->check_size();
        buffer_[this->size_-1] = c;
        return *this;
    }
    static_string& operator+=(const char* rhs)
    {
        const std::size_t last = this->size_;
        const std::size_t len  = traits_type::length(rhs);
        this->size_ += len;
        this->check_size();
        traits_type::copy(std::addressof(buffer_[last]), rhs, len+1);
        return *this;
    }
    template<std::size_t M>
    static_string& operator+=(const static_string<M>& rhs)
    {
        const std::size_t last = this->size_;
        this->size_ += rhs.size();
        this->check_size();
        traits_type::copy(std::addressof(buffer_[last]), rhs.c_str(), rhs.size()+1);
        return *this;
    }
    static_string& operator+=(const std::string& rhs)
    {
        const std::size_t last = this->size_;
        this->size_ += rhs.size();
        this->check_size();
        traits_type::copy(std::addressof(buffer_[last]), rhs.c_str(), rhs.size()+1);
        return *this;
    }

    bool empty() const noexcept {return this->size_ == 0;}
    void clear() {this->size_ = 0;}
    void resize(size_type i) {this->size_ = i; this->check_size();}

    value_type& operator[](size_type i)       noexcept {return buffer_[i];}
    value_type  operator[](size_type i) const noexcept {return buffer_[i];}
    value_type& at(size_type i)       {return buffer_.at(i);}
    value_type  at(size_type i) const {return buffer_.at(i);}

    value_type& front()       noexcept {return buffer_.front();}
    value_type  front() const noexcept {return buffer_.front();}
    value_type& back()        noexcept {return buffer_[this->size_-1];}
    value_type  back()  const noexcept {return buffer_[this->size_-1];}

    std::size_t size()     const noexcept {return this->size_;}
    std::size_t length()   const noexcept {return this->size_;}
    std::size_t capacity() const noexcept {return capacity_;}
    std::size_t max_size() const noexcept {return capacity_;}

    pointer       data()        noexcept {return buffer_.data();}
    const_pointer data()  const noexcept {return buffer_.data();}
    pointer       c_str()       noexcept {return buffer_.data();}
    const_pointer c_str() const noexcept {return buffer_.data();}

    iterator        begin()       noexcept {return buffer_.begin();}
    iterator        end()         noexcept {return buffer_.begin();}
    const_iterator  begin() const noexcept {return buffer_.begin();}
    const_iterator  end()   const noexcept {return buffer_.begin();}
    const_iterator cbegin() const noexcept {return buffer_.cbegin();}
    const_iterator cend()   const noexcept {return buffer_.cbegin();}

    reverse_iterator        rbegin()       noexcept {return buffer_.rbegin();}
    reverse_iterator        rend()         noexcept {return buffer_.rbegin();}
    const_reverse_iterator  rbegin() const noexcept {return buffer_.rbegin();}
    const_reverse_iterator  rend()   const noexcept {return buffer_.rbegin();}
    const_reverse_iterator crbegin() const noexcept {return buffer_.crbegin();}
    const_reverse_iterator crend()   const noexcept {return buffer_.crbegin();}

  private:

    std::size_t    size_;
    container_type buffer_;
};

template<std::size_t N>
constexpr std::size_t static_string<N>::capacity_;

template<std::size_t N, std::size_t M>
static_string<N+M>
operator+(const static_string<N>& lhs, const static_string<M>& rhs)
{
    static_string<N+M> retval(lhs);
    retval += rhs;
    return retval;
}
template<std::size_t N>
std::string operator+(const std::string& lhs, const static_string<N>& rhs)
{
    std::string retval(lhs);
    retval.append(rhs.c_str());
    return retval;
}
template<std::size_t N>
std::string operator+(const static_string<N>& lhs, const std::string& rhs)
{
    std::string retval(lhs.c_str());
    retval.append(rhs);
    return retval;
}
template<std::size_t N>
std::string operator+(const char* lhs, const static_string<N>& rhs)
{
    std::string retval(lhs);
    retval.append(rhs.c_str());
    return retval;
}
template<std::size_t N>
std::string operator+(const static_string<N>& lhs, const char* rhs)
{
    std::string retval(lhs.c_str());
    retval.append(rhs);
    return retval;
}

// append static_string into std::string
template<std::size_t N>
std::string& operator+=(std::string& lhs, const static_string<N>& rhs)
{
    lhs.append(rhs.c_str());
    return lhs;
}


template<std::size_t N, std::size_t M>
bool operator==(const static_string<N>& lhs, const static_string<M>& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) == 0;
}
template<std::size_t N, std::size_t M>
bool operator!=(const static_string<N>& lhs, const static_string<M>& rhs)
{
    return !(lhs == rhs);
}
template<std::size_t N, std::size_t M>
bool operator< (const static_string<N>& lhs, const static_string<M>& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) < 0;
}
template<std::size_t N, std::size_t M>
bool operator<=(const static_string<N>& lhs, const static_string<M>& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) <= 0;
}
template<std::size_t N, std::size_t M>
bool operator> (const static_string<N>& lhs, const static_string<M>& rhs)
{
    return !(lhs <= rhs);
}
template<std::size_t N, std::size_t M>
bool operator>=(const static_string<N>& lhs, const static_string<M>& rhs)
{
    return !(lhs < rhs);
}

template<std::size_t N>
std::ostream& operator<<(std::ostream& os, const static_string<N>& str)
{
    os << str.c_str();
    return os;
}

template<std::size_t N>
std::istream& operator>>(std::istream& is, static_string<N>& str)
{
    std::string tmp;
    is >> tmp;
    str = tmp;
    return is;
}

// compare with std::string.

template<std::size_t N>
bool operator==(const static_string<N>& lhs, const std::string& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) == 0;
}
template<std::size_t N>
bool operator==(const std::string& lhs, const static_string<N>& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) == 0;
}

template<std::size_t N>
bool operator!=(const static_string<N>& lhs, const std::string& rhs)
{
    return !(lhs == rhs);
}
template<std::size_t N>
bool operator!=(const std::string& lhs, const static_string<N>& rhs)
{
    return !(lhs == rhs);
}

template<std::size_t N>
bool operator< (const static_string<N>& lhs, const std::string& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) < 0;
}
template<std::size_t N>
bool operator< (const std::string& lhs, const static_string<N>& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) < 0;
}

template<std::size_t N>
bool operator<=(const static_string<N>& lhs, const std::string& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) <= 0;
}
template<std::size_t N>
bool operator<=(const std::string& lhs, const static_string<N>& rhs)
{
    return std::char_traits<char>::compare(lhs.data(), rhs.data(),
            std::min(lhs.size(), rhs.size())) <= 0;
}

template<std::size_t N>
bool operator> (const static_string<N>& lhs, const std::string& rhs)
{
    return !(lhs <= rhs);
}
template<std::size_t N>
bool operator> (const std::string& lhs, const static_string<N>& rhs)
{
    return !(lhs <= rhs);
}

template<std::size_t N>
bool operator>=(const std::string& lhs, const static_string<N>& rhs)
{
    return !(lhs < rhs);
}
template<std::size_t N>
bool operator>=(const static_string<N>& lhs, const std::string& rhs)
{
    return !(lhs < rhs);
}

} // mjolnir
#endif // MJOLNIR_UTIL_STATIC_STRING_HPP
