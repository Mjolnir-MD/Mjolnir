#ifndef MJOLNIR_UTIL_RANGE
#define MJOLNIR_UTIL_RANGE
#include <iterator>

namespace mjolnir
{

template<typename Iterator>
struct range
{
    typedef Iterator iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;

    range(iterator b, iterator e) noexcept
        : begin_(b), end_(e), sz_(std::distance(b, e))
    {}

    range() noexcept : sz_(0){}
    ~range() = default;
    range(const range&) = default;
    range(range&&)      = default;
    range& operator=(const range&) = default;
    range& operator=(range&&)      = default;

    iterator begin() const noexcept {return begin_;}
    iterator end()   const noexcept {return end_;}
    reverse_iterator rbegin() const noexcept {return reverse_iterator(begin_);}
    reverse_iterator rend()   const noexcept {return reverse_iterator(end_);}

    std::size_t size() const noexcept {return sz_;}

    iterator begin_;
    iterator end_;
    std::size_t sz_;
};

template<typename Iterator>
inline range<Iterator> make_range(Iterator b, Iterator e)
{
    return range<Iterator>(b, e);
}

} // mjolnir
#endif// MJOLNIR_UTIL_RANGE
