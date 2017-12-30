#ifndef JARNGREIPR_UTIL_SPLIT_STRING_ITERATOR
#define JARNGREIPR_UTIL_SPLIT_STRING_ITERATOR
#include <algorithm>
#include <iterator>
#include <string>

namespace mjolnir
{
namespace detail
{
template<typename stringT, bool is_delim_char>
struct split_string_iterator;
template<typename charT, typename traitsT, typename allocT, bool is_delim_char>
struct split_string_range;

template<typename stringT>
struct split_string_iterator<stringT, true>
{
    typedef stringT string_type;
    typedef string_type value_type;
    typedef value_type const* pointer;
    typedef value_type const& reference;
    typedef std::ptrdiff_t    difference_type;
    typedef std::bidirectional_iterator_tag iterator_category;
    typedef typename string_type::traits_type::char_type char_type;
    typedef split_string_iterator<string_type, true> self_type;

    split_string_iterator(const value_type& src, const char_type delim)
        noexcept
        : begin_(src.begin()), end_(std::find(src.begin(), src.end(), delim)),
          source_(src), delimiter_(delim)
    {}

    split_string_iterator(const value_type& src, const char_type delim,
                          const typename value_type::const_iterator bg,
                          const typename value_type::const_iterator ed)
        noexcept
        : begin_(bg), end_(ed), source_(src), delimiter_(delim)
    {}

    ~split_string_iterator() = default;
    split_string_iterator(const split_string_iterator&) = default;
    split_string_iterator(split_string_iterator&&) = default;
    split_string_iterator& operator=(const split_string_iterator&) = default;
    split_string_iterator& operator=(split_string_iterator&&) = default;

    self_type& operator++() noexcept
    {
        if(this->end_ == this->source_.end())
        {
            this->begin_ = this->end_;
        }
        else
        {
            this->begin_ = std::next(this->end_);
            this->end_ = std::find(this->begin_, this->source_.end(), this->delimiter_);
        }
        return *this;
    }
    self_type& operator--() noexcept
    {
        if(this->begin_ == this->source_.begin())
        {
            this->end_ = this->begin_;
        }
        else
        {
            this->end_ = std::prev(this->begin_);
            this->begin_ = std::find(std::reverse_iterator<
                    typename value_type::const_iterator>(this->end_),
                    this->source_.crend(), this->delimiter_).base();
        }
        return *this;
    }

    self_type operator++(int) noexcept {auto tmp(*this); ++(*this); return tmp;}
    self_type operator--(int) noexcept {auto tmp(*this); ++(*this); return tmp;}

    value_type operator*() const {return value_type(begin_, end_);}

    bool operator==(const split_string_iterator& rhs) const noexcept
    {
        return this->delimiter_ == rhs.delimiter_ && this->begin_ == rhs.begin_;
    }
    bool operator!=(const split_string_iterator& rhs) const noexcept
    {
        return !(*this == rhs);
    }

//   private:

    typename value_type::const_iterator begin_, end_;
    const value_type& source_;
    const char_type   delimiter_;
};

template<typename charT, typename traitsT, typename allocT>
struct split_string_range<charT, traitsT, allocT, true>
{
    typedef charT char_type;
    typedef std::basic_string<charT, traitsT, allocT> string_type;
    typedef split_string_iterator<string_type, true> iterator;
    typedef iterator const_iterator;
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    split_string_range(const string_type& source, const char_type delim)
        : src_(source), dlm_(delim)
    {}

    const_iterator begin()  const noexcept
    {return const_iterator(src_, dlm_);}
    const_iterator end()    const noexcept
    {return const_iterator(src_, dlm_, src_.end(), src_.end());}
    const_iterator cbegin() const noexcept
    {return const_iterator(src_, dlm_);}
    const_iterator cend()   const noexcept
    {return const_iterator(src_, dlm_, src_.end(), src_.end());}

  private:
    const string_type& src_;
    const char_type    dlm_;
};

} // detail

template<typename charT, typename traitsT = std::char_traits<charT>,
         typename allocT = std::allocator<charT>>
inline detail::split_string_range<charT, traitsT, allocT, true>
split_string(const std::basic_string<charT, traitsT, allocT>& str, charT delim)
{
    return detail::split_string_range<charT, traitsT, allocT, true>(str, delim);
}

} // mjolnir
#endif// JARNGREIPR_UTIL_SPLIT_STRING_ITERATOR
