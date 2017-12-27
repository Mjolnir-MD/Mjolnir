#include <jarngreipr/geninp/geninp.hpp>
#include <mjolnir/math/Vector.hpp>
#include <algorithm>
#include <map>

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        std::cerr << "Usage: $ jarngreipr <mode> [file.toml]\n";
        std::cerr << "mode geninp: generate input file\n";
        std::cerr << "             see input/example.toml for the format\n";
        std::cerr << std::flush;
        return 1;
    }

    const std::string mode(argv[1]);
    if(mode == "geninp")
    {
        return mjolnir::geninp<mjolnir::Vector<double, 3>>(--argc, ++argv);
    }
    else
    {
        std::cerr << "unknown mode: " << mode << '\n';
        std::cerr << "Usage: $ jarngreipr <mode> [file.toml]\n";
        std::cerr << "mode geninp: generate input file\n";
        std::cerr << "             see input/example.toml for the format\n";
        std::cerr << std::flush;
        return 1;
    }
    return 0;
}
