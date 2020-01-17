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
  public:

    progress_bar() noexcept
        : total_(1), r_total_(1.0)
    {}
    explicit progress_bar(std::size_t tot) noexcept
        : total_(tot), r_total_(1.0 / tot)
    {}
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

    void format(std::size_t count, std::ostream& os)
    {
        //XXX: requires UTF-8. TODO: consider setting locale
        constexpr auto full  = u8"█"; // U+2588 Full block
        constexpr auto l7    = u8"▉"; // U+2589 Left seven eighths block
        constexpr auto l6    = u8"▊"; // U+258A Left three quarters block
        constexpr auto l5    = u8"▋"; // U+258B Left five eighths block
        constexpr auto l4    = u8"▌"; // U+258C Left half block
        constexpr auto l3    = u8"▍"; // U+258D Left three eighths block
        constexpr auto l2    = u8"▎"; // U+258E Left one quarter block
        constexpr auto l1    = u8"▏"; // U+258F Left one eighth block

        const double ratio = (count == total_) ? 1.0 :
            std::max(0.0, std::min(1.0, count * this->r_total_));

        std::array<char, 8> buf1;
        std::snprintf(buf1.data(), 8, "%5.1f%%|", ratio * 100.0);
        os << '\r' << buf1.data();

        const std::size_t filled = std::floor(ratio*Width);
        for(std::size_t i=0; i<filled; ++i)
        {
            os << full;
        }
        if(Width > filled)
        {
            switch(static_cast<std::size_t>(ratio * Width * 8) % 8)
            {
                case 0:{os << ' '; break;}
                case 1:{os << l1;  break;}
                case 2:{os << l2;  break;}
                case 3:{os << l3;  break;}
                case 4:{os << l4;  break;}
                case 5:{os << l5;  break;}
                case 6:{os << l6;  break;}
                case 7:{os << l7;  break;}
            }
            for(std::size_t i=1; i<Width - filled; ++i)
            {
                os << ' ';
            }
        }
        os << '|' << std::flush;
        return ;
    }

    std::size_t total() const noexcept {return total_;}

  private:

    std::size_t total_;
    double      r_total_;
};

} // mjolnir
#endif// MJOLNIR_PROGRESS_BAR_HPP
