#include "base/vector.h"
#include "config/command_line.h"
#include "config/ArgManager.h"
#include "tools/Random.h"

EMP_BUILD_CONFIG( EvoConfig,
    GROUP(DEFAULT_GROUP, "General Settings"),
    VALUE(RANDOM_SEED, int, 0, "Random number seed (0 for based on time)"),
    VALUE(GENERATIONS, int, 10, "Number of generations to run for"),
    VALUE(N_GENOTYPES, int, 2, "Number of possible genotypes (must be >= 2)"),
    VALUE(K, double, 1000, "Carrying capacity"),
    VALUE(DEATH_RATE, double, .05, "Death rate (d)"),
    VALUE(MAX_BIRTH_RATE, double, 2, "Maximum birth rate (b0)")
)

int main(int argc, char* argv[])
{
    EvoConfig config;
    auto args = emp::cl::ArgManager(argc, argv);
    if (args.ProcessConfigOptions(config, std::cout, "MemicConfig.cfg", "Memic-macros.h") == false) exit(0);
    if (args.TestUnknown() == false) exit(0);  // If there are leftover args, throw an error.

    // Write to screen how the experiment is configured
    std::cout << "==============================" << std::endl;
    std::cout << "|    How am I configured?    |" << std::endl;
    std::cout << "==============================" << std::endl;
    config.Write(std::cout);
    std::cout << "==============================\n" << std::endl;

    emp::Random rnd(config.RANDOM_SEED());

    int N_GENOTYPES = config.N_GENOTYPES();

    // Set-up per-genotype values

    // Fitnesses (s in the equation) are relative to first genotype 
    // The first genotype doesn't need an s value
    emp::vector<double> rel_fitnesses;
    rel_fitnesses.resize(N_GENOTYPES - 1);

    // Initial population sizes of each genotype
    emp::vector<double> init_pops;
    init_pops.resize(N_GENOTYPES);

    // Outer vector is starting genotype, inner is to
    // e.g. mut_rates[3][4] is the probability of mutating
    // from genotype 3 to genotype 4.
    emp::vector<emp::vector<double> > mut_rates;
    mut_rates.resize(N_GENOTYPES);
    for (emp::vector<double> & vec : mut_rates) {
        vec.resize(N_GENOTYPES);
    }

    for (int i = 0; i < config.GENERATIONS(); i++) {
        
    }

}
 