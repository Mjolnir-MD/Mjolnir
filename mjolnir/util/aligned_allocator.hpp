#ifndef MJOLNIR_UTIL_ALIGNED_ALLOCATOR_HPP
#define MJOLNIR_UTIL_ALIGNED_ALLOCATOR_HPP
#include <type_traits>
#include <memory>
#include <limits>
#include <cstdlib>
#include <cassert>

#if defined(_WIN32)
#include <malloc.h>
#endif

namespace mjolnir
{

#if _POSIX_C_SOURCE >= 200112L

inline void* aligned_alloc(std::size_t alignment, std::size_t size)
{
    void *ptr = nullptr;
    if(posix_memalign(&ptr, alignment, size) != 0){return nullptr;}
    return ptr;
}

inline void aligned_free(void* ptr)
{
    std::free(ptr);
    return;
}

#elif defined(_WIN32)

inline void* aligned_alloc(std::size_t alignment, std::size_t size)
{
    return _aligned_malloc(size, alignment);
}

inline void aligned_free(void* ptr)
{
    _aligned_free(ptr);
    return;
}

#else // fallback...

// This implementation uses C++11 std::align() that finds an aligned region from
// larger space. Since it changes a pointer to the start of the allocated area,
// we need to store the original pointer to correctly free the region. To store
// the original pointer, it allocates an extra region for 1x `void*`. Because
// `aligned_free` does not know the length of the region, we cannot put the
// original pointer at the end of the region. The only place we can put is
// before the starting point of the aligned region. Since the size of padding
// is also unknown `from aligned_free`, we need to put the value just before
// the aligned region (of which pointer would be passed to `aligned_free`).
// Thus the resulting memory region looks like the following.
//
//  ______________ allocated region ____________________
// '                                                    '
// |(unknown size)|sizeof(pointer) | required size      |
// |(padding)     |original pointer| aligned region ... |
// ^               |               ^- pointer to this location will be returned.
// +---------------+ this pointer (returned from malloc) points the area itself.
//
// So we need to allocate a memory region that has size
// `required size + alignment padding + sizeof(void*)`.

inline void* aligned_alloc(std::size_t alignment, std::size_t size)
{
    // malloc returns a memory region that is aligned for any scalar type.
    // So at least the alignment of void* is guaranteed.
    constexpr std::size_t minimum_align = alignof(void*);

    if(alignment <= minimum_align)
    {
        // alignment is automatically guaranteed.
        // for the consistency (for aligned_free), it stores the original ptr.
        void* ptr = std::malloc(size + sizeof(void*));
        *(reinterpret_cast<void**>(ptr)) = ptr;
        return reinterpret_cast<void*>(reinterpret_cast<void**>(ptr) + 1);
    }
    else
    {
        assert(alignment > minimum_align);
        const std::size_t offset = alignment - minimum_align;
        void* const ptr = std::malloc(size + sizeof(void*) + offset);

        // prepare space to write the original ptr.
        void* aligned_ptr =
            reinterpret_cast<void*>(reinterpret_cast<void**>(ptr) + 1);
                
        // search aligned region ...
        std::size_t space = size + offset;
        void* tmp = std::align(alignment, size, aligned_ptr, space);
        assert(tmp == aligned_ptr);

        // write the original pointer
        *(reinterpret_cast<void**>(aligned_ptr) - 1) = ptr;
        return aligned_ptr;
    }
}

inline void aligned_free(void* ptr)
{
    // free the region that is pointed by the original (before-aligned) pointer
    if(ptr)
    {
        std::free(*(reinterpret_cast<void**>(ptr) - 1));
    }
    return;
}
#endif // aligned_alloc/free

template<typename T, std::size_t Alignment = std::alignment_of<T>::value>
class aligned_allocator
{
  public:
    using value_type      = T;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;
    using pointer         = value_type*;
    using const_pointer   = value_type const*;
    using reference       = value_type&;
    using const_reference = value_type const&;
    using propagate_on_container_move_assignment = std::true_type;

    // use the maximum alignment in {Alignment, alignof(T), alignof(void*)}.
    static constexpr std::size_t alignment = compiletime::max(
            Alignment, compiletime::max(alignof(T), alignof(void*))
        );

    template<typename U>
    struct rebind
    {
        using other = aligned_allocator<U, alignment>;
    };

  public:

    aligned_allocator() noexcept = default;
    aligned_allocator(const aligned_allocator&) noexcept = default;
    aligned_allocator& operator=(const aligned_allocator&) noexcept = default;

    template<typename U, std::size_t B>
    explicit aligned_allocator(const aligned_allocator<U, B>&) noexcept {}
    template<typename U, std::size_t B>
    aligned_allocator& operator=(const aligned_allocator<U, B>&) noexcept
    {return *this;}

    pointer allocate(std::size_t n)
    {
        void* ptr = aligned_alloc(alignment, sizeof(T) * n);
        if(!ptr) {throw std::bad_alloc{};}
        return reinterpret_cast<pointer>(ptr);
    }
    void deallocate(pointer p, std::size_t)
    {
        aligned_free(p);
    }

    pointer       address(reference       x) const noexcept
    {return std::addressof(x);}
    const_pointer address(const_reference x) const noexcept
    {return std::addressof(x);}

    size_type max_size() const noexcept
    {
        return std::numeric_limits<size_type>::max() /
               std::max(sizeof(value_type), alignment);
    }

    void construct(pointer p, const_reference val)
    {
        new(reinterpret_cast<void*>(p)) value_type(val);
        return;
    }
    template <class U, class... Args>
    void construct(U* p, Args&&... args)
    {
        new(reinterpret_cast<void*>(p)) U(std::forward<Args>(args)...);
        return;
    }

    void destroy(pointer p){p->~value_type(); return;}
    template <class U>
    void destroy(U* p){p->~U(); return;}

};
template<typename T, std::size_t A>
constexpr std::size_t aligned_allocator<T, A>::alignment;

template<typename T, std::size_t A>
inline bool
operator==(const aligned_allocator<T, A>&, const aligned_allocator<T, A>&)
{
    return true;
}
template<typename T, std::size_t A>
inline bool
operator!=(const aligned_allocator<T, A>&, const aligned_allocator<T, A>&)
{
    return false;
}

} // mjolnir
#endif // MJOLNIR_ALIGNED_ALLOCATOR_HPP
