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

    while(simulator->step()){/* do nothing */;}

    simulator->finalize();

    const auto stop = std::chrono::system_clock::now();
    std::cout << "elapsed time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                      stop - start).count() << " [msec]" << std::endl;

    return 0;
}
