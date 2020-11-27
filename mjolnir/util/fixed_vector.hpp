#ifndef MJOLNIR_UTIL_FIXED_VECTOR_HPP
#define MJOLNIR_UTIL_FIXED_VECTOR_HPP
#include <algorithm>
#include <array>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <string>
#include <mjolnir/util/string.hpp>

namespace mjolnir
{

template<typename T, std::size_t N>
struct fixed_vector
{
  public:
    using value_type             = T;
    using reference              = value_type&;
    using const_reference        = value_type const&;
    using pointer                = value_type*;
    using const_pointer          = value_type const*;
    using size_type              = std::size_t;
    using difference_type        = std::ptrdiff_t;
    using backend_type           = std::array<T, N>;
    using iterator               = typename backend_type::iterator;
    using const_iterator         = typename backend_type::const_iterator;
    using reverse_iterator       = typename backend_type::reverse_iterator;
    using const_reverse_iterator = typename backend_type::const_reverse_iterator;

// XXX This requires T to be default constructible. To remove this requirement,
//     we need to replace std::array by a correctly aligned fixed-size storage
//     and use placement new to construct a value on it. But currently we only
//     use this struct for default constructible types, we keep it simple.
//     ... This also means that the default constructor of value_type always
//     runs N times at the construction, but we currently assume that the cost
//     is almost ignorable.
    static_assert(std::is_default_constructible<value_type>::value, "");

  public:

    fixed_vector() noexcept: size_(0) {}
    ~fixed_vector() = default;
    fixed_vector(fixed_vector const&) = default;
    fixed_vector(fixed_vector &&)     = default;
    fixed_vector& operator=(fixed_vector const&) = default;
    fixed_vector& operator=(fixed_vector &&)     = default;

    fixed_vector(const std::size_t sz) : size_(sz)
    {
        if(N < sz)
        {
            throw std::bad_alloc("fixed_vector(size_t): too many elements");
        }
    }
    fixed_vector(const std::size_t sz, const value_type& x): size_(sz)
    {
        if(N < sz)
        {
            throw std::bad_alloc("fixed_vector(size_t, value_type): too many elements");
        }
        std::fill_n(container_.begin(), sz, x);
    }
    template<typename InputIterator>
    fixed_vector(InputIterator first, InputIterator last)
    {
        const auto sz = std::distance(first, last);
        if(N < sz)
        {
            throw std::bad_alloc("fixed_vector(first, last): too many elements");
        }
        this->size_ = sz;
        std::copy(first, last, container_.begin());
    }
    fixed_vector(std::initializer_list<T> il)
        : size_(il.size())
    {
        if(N < il.size())
        {
            size_ = 0;
            throw std::bad_alloc("fixed_vector(initializer_list): too many elements");
        }
        std::copy(il.begin(), il.end(), container_.begin());
    }
    fixed_vector& operator=(std::initializer_list<T> il)
    {
        if(N < il.size())
        {
            throw std::bad_alloc("fixed_vector& operator=(initializer_list): too many elements");
        }
        size_ = il.size();
        std::copy(il.begin(), il.end(), container_.begin());
    }

    void push_back(const value_type& v)
    {
        if(this->size_ == N)
        {
            throw std::bad_alloc("fixed_vector::push_back(): too many elements");
        }
        container_[size_] = v;
    }
    void push_back(value_type&& v)
    {
        if(this->size_ == N)
        {
            throw std::bad_alloc("fixed_vector::push_back(): too many elements");
        }
        container_[size_] = std::move(v);
    }
    template<typename ... Ts>
    void emplace_back(Ts&& ... args)
    {
        if(this->size_ == N)
        {
            throw std::bad_alloc("fixed_vector::emplace_back(): too many elements");
        }
        container_[size_] = value_type(std::forward<Ts>(args)...);
    }
    void pop_back()
    {
        if(size_ != 0)
        {
            size_ -= 1;
        }
    }

    // TODO insert, emplace

    iterator erase(const_iterator pos)
    {
        if(pos == container_.end())
        {
            return container_.end();
        }
        assert(this->size_ != 0);

        this->size_ -= 1;
        std::copy(std::next(pos), container_.end(), pos);
    }
    iterator erase(const_iterator first, const_iterator last)
    {
        this->size_ -= std::distance(first, last);
        std::copy(last, container_.end(), first);
    }

    void swap(fixed_vector& other)
    {
        using std::swap;
        swap(other.size_,      this->size_);
        swap(other.container_, this->container_);
    }
    void clear() {size_ = 0;}

    value_type&       operator[](size_type i)       noexcept {return container_[i];}
    value_type const& operator[](size_type i) const noexcept {return container_[i];}
    value_type&       at(size_type i)
    {
        if(size_ <= i)
        {
            using namespace ::mjolnir::literals::string_literals;
            std::out_of_range("fixed_vector::at("_s + std::to_string(i) + "): "
                "index exceeds the current size (" + std::to_string(size_) + ")");
        }
        return container_.at(i);
    }
    value_type const& at(size_type i) const
    {
        if(size_ <= i)
        {
            using namespace ::mjolnir::literals::string_literals;
            std::out_of_range("fixed_vector::at("_s + std::to_string(i) + "): "
                "index exceeds the current size (" + std::to_string(size_) + ")");
        }
        return container_.at(i);
    }
    value_type*       data()        noexcept {return container_.data();}
    value_type const* data()  const noexcept {return container_.data();}
    value_type&       front()       noexcept {return container_.front();}
    value_type const& front() const noexcept {return container_.front();}
    value_type&       back()        noexcept {return container_[size_-1];}
    value_type const& back()  const noexcept {return container_[size_-1];}

    void resize(size_type sz)
    {
        if(N < sz)
        {
            throw std::bad_alloc("fixed_vector::resize()");
        }
        size_ = sz;
    }
    void resize(size_type sz, const value_type& c)
    {
        if(N < sz)
        {
            throw std::bad_alloc("fixed_vector::resize()");
        }
        if(size_ < sz)
        {
            std::fill_n(container_.begin() + size_, sz - size_, c);
        }
        size_ = sz;
    }

    std::size_t capacity() const noexcept {return N;}
    std::size_t max_size() const noexcept {return N;}
    std::size_t size()     const noexcept {return size_;}
    bool        empty()    const noexcept {return size_ == 0;}

    iterator       begin()        noexcept {return container_.begin();}
    iterator       end()          noexcept {return container_.begin() + size_;}
    const_iterator begin()  const noexcept {return container_.begin();}
    const_iterator end()    const noexcept {return container_.begin() + size_;}
    const_iterator cbegin() const noexcept {return container_.begin();}
    const_iterator cend()   const noexcept {return container_.begin() + size_;}

    reverse_iterator       rbegin()        noexcept {return reverse_iterator(this->end());   }
    reverse_iterator       rend()          noexcept {return reverse_iterator(this->begin()); }
    const_reverse_iterator rbegin()  const noexcept {return const_reverse_iterator(this->cend());  }
    const_reverse_iterator rend()    const noexcept {return const_reverse_iterator(this->cbegin());}
    const_reverse_iterator crbegin() const noexcept {return const_reverse_iterator(this->cend());  }
    const_reverse_iterator crend()   const noexcept {return const_reverse_iterator(this->cbegin());}

  private:
    std::size_t size_;
    std::array<T, N> container_;
};

} // mjolnir
#endif// MJOLNIR_UTIL_FIXED_VECTOR_HPP
