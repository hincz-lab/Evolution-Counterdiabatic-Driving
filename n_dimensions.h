#include "base/vector.h"
#include "config/command_line.h"
#include "config/ArgManager.h"
#include "tools/Random.h"
#include "tools/vector_utils.h"
#include "tools/IndexMap.h"

EMP_BUILD_CONFIG( EvoConfig,
    GROUP(DEFAULT_GROUP, "General Settings"),
    VALUE(RANDOM_SEED, int, 0, "Random number seed (0 for based on time)"),
    VALUE(GENERATIONS, int, 10, "Number of generations to run for"),
    VALUE(N_GENOTYPES, int, 2, "Number of possible genotypes (must be >= 2)"),
    VALUE(K, double, 1000, "Carrying capacity"),
    VALUE(DEATH_RATE, double, .05, "Death rate (d)"),
    VALUE(MAX_BIRTH_RATE, double, 2, "Maximum birth rate (b0)")
)

/*Heaviside Theta Function*/ 
int H(double z){
    return z>0;
}

class NDimSim {
    private:
    emp::Random rnd;

    // Vectors to hold per-genotype values
    emp::vector<double> rel_fitnesses;
    emp::vector<double> init_pops;
    emp::vector<double> current_pops;
    emp::vector<double> new_pops;
    emp::vector<emp::IndexMap> mut_rates;

    // Localized config parameters
    int N_GENOTYPES;
    int GENERATIONS;
    double K;
    double DEATH_RATE;
    double MAX_BIRTH_RATE;

    public:
    NDimSim(EvoConfig & config) {
        Setup(config);
    }

    void Setup(EvoConfig & config) {
        rnd = emp::Random(config.RANDOM_SEED());

        // Localize config parameters
        N_GENOTYPES = config.N_GENOTYPES();
        GENERATIONS = config.GENERATIONS();
        K = config.K();
        DEATH_RATE = config.DEATH_RATE();
        MAX_BIRTH_RATE = config.MAX_BIRTH_RATE();

        // Set-up per-genotype values

        // Fitnesses (s in the equation) are relative to first genotype 
        // The first genotype's s will be 0
        rel_fitnesses.resize(N_GENOTYPES);

        // Initial population sizes of each genotype
        init_pops.resize(N_GENOTYPES);

        // Current population sizes of each genotype
        current_pops.resize(N_GENOTYPES);

        // Holder for population sizes currently being calculated
        new_pops.resize(N_GENOTYPES);

        // Outer vector is starting genotype, inner index map is to
        // e.g. the bin size of mut_rates[3][4] represents the probability of mutating
        // from genotype 3 to genotype 4.
        // Diagonal is probability of keeping same genotype
        // An IndexMap is a data structure that effectively maintains
        // a set of "bins" of different weights, corresponding to a probability
        // of a given outcome. They are useful for efficiently stochastically selecting
        // a single outcome when the probabilities of all possibilites are known
        mut_rates.resize(N_GENOTYPES);
        for (emp::IndexMap vec : mut_rates) {
            vec.Resize(N_GENOTYPES);
        }

    }

    /* Birth rate of first genotype */
    /* Adapted from code written by Julia Pelesko */
    double Birth(){

        double T = emp::Sum(current_pops)/K;
        double m; 

        if (T > 1) {
            m = 0; 
        } else {
            m = MAX_BIRTH_RATE*(1-T); 
        }

        return m; 
    }

    /* Birth rate of other genotypes (get's called for 
    all genotypes, but will not modify the value returned
    by Birth() for the focal genotype because its relative
    fitness will be 0) */
    /* Adapted from code written by Julia Pelesko */
    double Birth(int genotype) {
        return Birth()/(1 + rel_fitnesses[genotype]);
    }

    void RunStep() {
        for (size_t i = 0; i < new_pops.size(); i++) {
            new_pops[i] = 0;
        }

        for (int genotype = 0; genotype < N_GENOTYPES; genotype++) {
            // Store reference to current probability map, for simplicity
            emp::IndexMap & mut_probs = mut_rates[genotype];
            emp_assert(mut_probs.GetWeight() == 1 && "Mutation probabilities should sum to 1");
            
            // Loop through all individuals of this genotype and figure out what happens to them
            for (int individual = 0; individual < current_pops[genotype]; individual++) {
                
                // Check if individual died
                if (rnd.P(1 - DEATH_RATE)) {
                    new_pops[genotype]++; // Individual is still alive so we count it

                    // Check if indiviudal gives birth
                    if (rnd.P(Birth(genotype))) {
                        // Select genotype of offsping based on transition probabilities
                        size_t new_genotype = mut_probs.Index(rnd.GetDouble(1));
                        new_pops[new_genotype]++;             
                    }
                }
            }
        }

    }

    void Run() {
        for (int gen = 0; gen < GENERATIONS; gen++) {
            RunStep();
        }        
    }
};
