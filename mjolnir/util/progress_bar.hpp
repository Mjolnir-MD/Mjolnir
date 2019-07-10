#ifndef MJOLNIR_UTIL_PROGRESS_BAR_HPP
#define MJOLNIR_UTIL_PROGRESS_BAR_HPP
#include <array>
#include <cstring>
#include <cmath>
#include <cassert>

namespace mjolnir
{

template<std::size_t Width>
class progress_bar
{
    //XXX: requires UTF-8. TODO: consider setting locale
    static constexpr auto full  = u8"█"; // U+2588 Full block
    static constexpr auto l7    = u8"▉"; // U+2589 Left seven eighths block
    static constexpr auto l6    = u8"▊"; // U+258A Left three quarters block
    static constexpr auto l5    = u8"▋"; // U+258B Left five eighths block
    static constexpr auto l4    = u8"▌"; // U+258C Left half block
    static constexpr auto l3    = u8"▍"; // U+258D Left three eighths block
    static constexpr auto l2    = u8"▎"; // U+258E Left one quarter block
    static constexpr auto l1    = u8"▏"; // U+258F Left one eighth block

    // \r100.0%| ... W*3 ... |\0
    static constexpr std::size_t buffer_size = 10 + Width * 3;

  public:

    progress_bar(): total_(1), r_total_(1.0)
    {
        buffer_.fill('\0');
        buffer_.front() = '\r';
    }
    explicit progress_bar(std::size_t tot) : total_(tot), r_total_(1.0 / tot)
    {
        buffer_.fill('\0');
        buffer_.front() = '\r';
    }
    ~progress_bar() = default;
    progress_bar(progress_bar const&) = default;
    progress_bar(progress_bar &&)     = default;
    progress_bar& operator=(progress_bar const&) = default;
    progress_bar& operator=(progress_bar &&)     = default;

    void reset(const std::size_t total_step)
    {
        this->total_   = total_step;
        this->r_total_ = 1.0 / static_cast<double>(total_step);
        return;
    }

    const char* format(std::size_t count)
    {
        const double ratio = (count == total_) ? 1.0 : count * this->r_total_;
        assert(ratio <= 1.0);

        char* iter = buffer_.data();
        iter++; // the first character is always \r
        iter += std::sprintf(iter, "%5.1f", ratio * 100.0);
        *iter++ = '%';
        *iter++ = '|';

        const std::size_t filled = std::floor(ratio*Width);
        for(std::size_t i=0; i<filled; ++i)
        {
            *iter++ = full[0];
            *iter++ = full[1];
            *iter++ = full[2];
        }
        if(Width > filled)
        {
            switch(static_cast<std::size_t>(ratio * Width * 8) % 8)
            {
                case 0:{*iter++ = ' '; break;}
                case 1:{*iter++ = l1[0]; *iter++ = l1[1]; *iter++ = l1[2]; break;}
                case 2:{*iter++ = l2[0]; *iter++ = l2[1]; *iter++ = l2[2]; break;}
                case 3:{*iter++ = l3[0]; *iter++ = l3[1]; *iter++ = l3[2]; break;}
                case 4:{*iter++ = l4[0]; *iter++ = l4[1]; *iter++ = l4[2]; break;}
                case 5:{*iter++ = l5[0]; *iter++ = l5[1]; *iter++ = l5[2]; break;}
                case 6:{*iter++ = l6[0]; *iter++ = l6[1]; *iter++ = l6[2]; break;}
                case 7:{*iter++ = l7[0]; *iter++ = l7[1]; *iter++ = l7[2]; break;}
            }
            for(std::size_t i=0; i<Width - filled - 1; ++i)
            {
                *iter++ = ' ';
            }
        }
        *iter++ = '|';
        *iter++ = '\0';
        return buffer_.data();
    }

    std::size_t total() const noexcept {return total_;}

  private:

    std::size_t total_;
    double      r_total_;
    std::array<char, buffer_size> buffer_;
};


} // mjolnir
#endif// MJOLNIR_PROGRESS_BAR_HPP
