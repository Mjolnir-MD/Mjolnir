#include <mjolnir/input/read_input_file.hpp>

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input.toml>" << std::endl;
        return 1;
    }

    std::cerr << "reading input file...\n";
    auto simulator = mjolnir::read_input_file(std::string(argv[1]));
    std::cerr << "done." << std::endl;

    std::cerr << "initializing simulator...\n";
    simulator->initialize();
    std::cerr << "done." << std::endl;

    const auto start = std::chrono::system_clock::now();
    std::cerr << "start running simulation" << std::endl;

    simulator->run();

    simulator->finalize();

    const auto stop = std::chrono::system_clock::now();
    const auto total = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::cerr << "elapsed time: ";
    std::cerr << std::fixed << std::setprecision(1);
    if(total < 1000.0)
    {
        std::cerr << total << " [msec]";
    }
    else if(total < 1000.0 * 60.0)
    {
        std::cerr << total * 0.001 << " [sec]";
    }
    else if(total < 1000.0 * 60.0 * 60.0)
    {
        std::cerr << total * 0.001 * 0.0167 << " [min]";
    }
    else if(total < 1000.0 * 60.0 * 60.0 * 24.0)
    {
        std::cerr << total * 0.001 * 0.0167 * 0.0167 << " [hr]";
    }
    else
    {
        std::cerr << total * 0.001 * 0.0167 * 0.0167 * 0.0417 << " [day]";
    }
    std::cerr << std::endl;

    return 0;
}
