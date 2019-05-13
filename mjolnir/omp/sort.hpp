#ifndef MJOLNIR_OMP_SORT_HPP
#define MJOLNIR_OMP_SORT_HPP
#include <type_traits>
#include <algorithm>
#include <functional>
#include <utility>
#include <vector>
#include <omp.h>

namespace mjolnir
{
namespace omp
{

// Do NOT call this from parallel region.
template<typename Value, typename Alloc, typename Comparator>
void sort(std::vector<Value, Alloc>& vec, std::vector<Value, Alloc>& buf,
          Comparator comp)
{
    const std::size_t max_threads = omp_get_max_threads();
    if(max_threads == 1 || vec.size() < max_threads * 4)
    {
        std::sort(vec.begin(), vec.end(), comp);
        return;
    }
    if(vec.size() != buf.size())
    {
        buf.resize(vec.size());
    }

#pragma omp parallel shared(vec, buf, comp)
    {
        const auto num_threads   = omp_get_num_threads();
        const auto thread_id     = omp_get_thread_num();
        const auto range_size    = vec.size();
        const auto subrange_size = (range_size % num_threads == 0) ?
            (range_size / num_threads) : (range_size / num_threads + 1);
        const auto offset = thread_id * subrange_size;

        if(offset < range_size)
        {
            const auto first = vec.begin() + offset;
            const auto last  = vec.begin() +
                std::min(offset + subrange_size, range_size);
            std::sort(first, last, comp);
        }
#pragma omp barrier

        int num_ranges = num_threads;
        int delta      = 2;
        int first      = offset;
        int middle     = offset + subrange_size;
        int last       = std::min(offset + 2 * subrange_size, range_size);
        while(num_ranges > 1)
        {
            if(thread_id % delta == 0 && middle < last)
            {
                std::merge(vec.begin() + first,  vec.begin() + middle,
                           vec.begin() + middle, vec.begin() + last,
                           buf.begin() + first,  comp);

                //     merged        merged (by another thread)
                // .------+------.------+------.
                // |......|......|......|......|
                // f      m      l -----+
                // f      +-->  m2      +----> l2
                middle  = last;
                last   += delta * subrange_size;
                last    = std::min(last, static_cast<int>(range_size));
            }
            else if(thread_id % delta == 0)
            {
                std::copy(vec.begin() + first, vec.begin() + last, buf.begin() + first);
            }

            delta *= 2;
            num_ranges = (num_ranges + 1) / 2;
#pragma omp barrier
#pragma omp single
            {
                // both are std::vector. we don't need to use ADL.
                std::swap(vec, buf);
            } // an implicit barrier here
        }
    }
    return;
}

template<typename Value, typename Alloc>
void sort(std::vector<Value, Alloc>& vec, std::vector<Value, Alloc>& buf)
{
    using comparator_type = std::less<Value>;
    ::mjolnir::omp::sort(vec, buf, comparator_type());
    return;
}

} // omp
} // mjolnir
#endif // MJOLNIR_OMP_SORT_HPP
