#ifndef MJOLNIR_UTIL_RANGE_HPP
#define MJOLNIR_UTIL_RANGE_HPP
#include <mjolnir/util/throw_exception.hpp>
#include <iterator>

namespace mjolnir
{

template<typename Iterator>
struct range
{
    using iterator         = Iterator;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using value_type       = typename std::iterator_traits<iterator>::value_type;
    using reference        = typename std::iterator_traits<iterator>::reference;

    static_assert(std::is_same<std::random_access_iterator_tag,
        typename std::iterator_traits<iterator>::iterator_category
        >::value, "mjolnir::range assumes iterator to be randomly accessible");

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

    reference operator[](const std::size_t i) const noexcept
    {
        return *(this->begin_ + i);
    }
    reference at(const std::size_t i) const noexcept
    {
        if(i >= sz_)
        {
            throw_exception<std::out_of_range>("mjolnir::range::at(i = ", i,
                    "), index out of range (this->size() = ", this->sz_, ")");
        }
        return *(this->begin_ + i);
    }

    std::size_t size() const noexcept {return sz_;}
    bool empty() const noexcept {return sz_ == 0u;}

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
#endif// MJOLNIR_UTIL_RANGE_HPP
