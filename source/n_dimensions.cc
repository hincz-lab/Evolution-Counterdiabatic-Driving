#include "n_dimensions.h"

int main(int argc, char* argv[])
{
    EvoConfig config;
    config.Read("NDim.cfg");
    auto args = emp::cl::ArgManager(argc, argv);
    if (args.ProcessConfigOptions(config, std::cout, "NDim.cfg", "NDim-macros.h") == false) exit(0);
    if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

    // Write to screen how the experiment is configured
    std::cout << "==============================" << std::endl;
    std::cout << "|    How am I configured?    |" << std::endl;
    std::cout << "==============================" << std::endl;
    config.Write(std::cout);
    std::cout << "==============================\n" << std::endl;

    NDimSim sim(config);
    sim.Run();

}
 