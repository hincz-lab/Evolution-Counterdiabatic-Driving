#include "base/vector.h"
#include "config/command_line.h"
#include "config/ArgManager.h"
#include "tools/Random.h"
#include "tools/vector_utils.h"
#include "tools/IndexMap.h"
#include "tools/File.h"
#include "data/DataFile.h"

EMP_BUILD_CONFIG( EvoConfig,
    GROUP(DEFAULT_GROUP, "General Settings"),
    VALUE(RANDOM_SEED, int, 0, "Random number seed (0 for based on time)"),
    VALUE(GENERATIONS, int, 10, "Number of generations to run for"),
    VALUE(N_GENOTYPES, int, 2, "Number of possible genotypes (must be >= 2)"),
    VALUE(K, double, 1000, "Carrying capacity"),
    VALUE(DEATH_RATE, double, .05, "Death rate (d)"),
    VALUE(MAX_BIRTH_RATE, double, 2, "Maximum birth rate (b0)"),

    GROUP(PER_GENOTYPE_VALUES, "Per-genotype values"),
    VALUE(FITNESSES, std::string, "0,1", "Either a list of relative fitnesses, separated by commas, or a file containing them"),
    VALUE(INIT_POPS, std::string, "100,10", "Either a list of initial population sizes, separated by commas, or a file containing them"),
    VALUE(TRANSITION_PROBS, std::string, ".95,.05:.05,.95", 
        "Either a matrix of transition probabilities or a file containing one. Rows are original genotype, columns are new one. Use commas to separate values within rows. In files, use newlines between rows. On command-line, use colons.")
)

/*Heaviside Theta Function*/ 
int H(double z){
    return z>0;
}

// Functions for extracting parameter lists from files
emp::File ExtractStringFromFile(std::string & fname, std::string param_name) {
    std::cout << "Loading " << param_name << " from file: " << fname <<
        ". If you were not trying to load a file, check for stray " << 
        "letters in your " << param_name << " list." << std::endl;
    emp::File param_file(fname);
    if (param_file.GetNumLines() > 1) {
        std::cout << 
            "Warning: file containing " << param_name << " has more than one "
            << "line. The expected format is a single line with "
            << param_name << " specified as numbers and separated by " 
            << "commas." << std::endl;
    }

    return param_file;
}

std::string ExtractStringFromOneLineFile(std::string & fname, std::string param_name) {
    emp::File param_file = ExtractStringFromFile(fname, param_name);
    return param_file.front();
}

emp::vector<std::string> ExtractMultiLineFile(std::string & fname, std::string param_name) {
    emp::File param_file = ExtractStringFromFile(fname, param_name);
    return param_file.GetAllLines();
}


class NDimSim {
    private:
    emp::Random rnd; // Random-number generator
    emp::DataFile pop_sizes; // For tracking pops over time
    int curr_gen = 0; // Current generation

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
    std::string FITNESSES;
    std::string INIT_POPS;
    std::string TRANSITION_PROBS;

    public:
    NDimSim(EvoConfig & config) : pop_sizes("pop_sizes.csv") {
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
        FITNESSES = config.FITNESSES();
        INIT_POPS = config.INIT_POPS();
        TRANSITION_PROBS = config.TRANSITION_PROBS();

        // Set-up per-genotype values

        // Fitnesses (s in the equation) are relative to first genotype 
        // The first genotype's s will be 0
        rel_fitnesses.resize(N_GENOTYPES);
        InitializeFitnesses();

        // Initial population sizes of each genotype
        init_pops.resize(N_GENOTYPES);
        InitializeInitPops();

        // Current population sizes of each genotype
        current_pops = init_pops;

        // Holder for population sizes currently being calculated
        // Gets zeroed out in RunStep();
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
        for (emp::IndexMap & vec : mut_rates) {
            vec.Resize(N_GENOTYPES);
        }
        InitializeTransitionProbs();

        // Set up data tracking
        pop_sizes.AddVar(curr_gen,"generation");
        for (int i = 0; i < N_GENOTYPES; i++) {
            // Using a reference to vector contents is safe here
            // because vector size will never change
            pop_sizes.AddVar(current_pops[i],"pop" + emp::to_string(i));
        }
        pop_sizes.PrintHeaderKeys();
        pop_sizes.SetTimingRepeat(1);


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
                    if (Birth(genotype) - rnd.GetDouble(1) >= 0) {
                        // Select genotype of offsping based on transition probabilities
                        size_t new_genotype = mut_probs.Index(rnd.GetDouble(1));
                        new_pops[new_genotype]++;             
                    }
                }
            }
        }

        // Update current population sizes
        std::swap(current_pops, new_pops);

    }

    void Run() {
        for (int gen = 0; gen < GENERATIONS; gen++) {
            curr_gen = gen;
            pop_sizes.Update(gen);
            RunStep();
        }
        curr_gen = GENERATIONS; // So that time stamp in last line is correct
        pop_sizes.Update(GENERATIONS); // Record final counts
    }

    void InitializeFitnesses() {
        if (emp::has_letter(FITNESSES)) {
            // FITNESSES is a file
            FITNESSES = ExtractStringFromOneLineFile(FITNESSES, "fitness");
        }

        emp::vector<std::string> sliced_fits = emp::slice(FITNESSES, ',');

        if ((int)sliced_fits.size() != N_GENOTYPES) {
            if ((int)sliced_fits.size() == N_GENOTYPES - 1) {
                // Assume missing fitness is the one for the
                // genotype all other fitnesses are relative to
                // (the first genotype)
                sliced_fits.insert(sliced_fits.begin(),"0"); 
            } else {
                std::cout << "Error: Not enough fitnesses supplied." << 
                    " Attempting to assign " << sliced_fits.size() << 
                    " to " << N_GENOTYPES << " genotypes." << std::endl;
                exit(1);
            }
        }

        if (sliced_fits[0] != "0") {
            std::cout << "Warning: Relative fitness of focal genotype is not 0" << std::endl;
        }

        std::cout << "Relative fitnesses: " << std::endl;
        for (int i = 0; i < N_GENOTYPES; i++) {
            rel_fitnesses[i] = emp::from_string<double>(sliced_fits[i]);
            std::cout << i << ": " << rel_fitnesses[i] << std::endl;
        }
        std::cout << std::endl;
    }

    void InitializeInitPops() {
        if (emp::has_letter(INIT_POPS)) {
            // INIT_POPS is a filename
            INIT_POPS = ExtractStringFromOneLineFile(INIT_POPS, "initial population");
        }

        emp::vector<std::string> sliced_pops = emp::slice(INIT_POPS, ',');

        if ((int)sliced_pops.size() != N_GENOTYPES) {
            std::cout << "Error: Not enough initial population sizes supplied." << 
                " Attempting to assign " << sliced_pops.size() << 
                " to " << N_GENOTYPES << " genotypes." << std::endl;
            exit(1);
        }

        std::cout << "Initial populations: " << std::endl;
        for (int i = 0; i < N_GENOTYPES; i++) {
            init_pops[i] = emp::from_string<double>(sliced_pops[i]);
            std::cout << i << ": " << init_pops[i] << std::endl;
        }
        std::cout << std::endl;
    }

    void InitializeTransitionProbs() {
        emp::vector<std::string> rows;
        if (emp::has_one_of(TRANSITION_PROBS, "abcdefghijklmopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")) {
            // INIT_POPS is a filename
            // This will technically fail if the file name is composed excelusively of 
            // lower case ns, numbers, and punctuation, but that seems unlikely.
            rows = ExtractMultiLineFile(TRANSITION_PROBS, "transition probabilities");
        } else {
            rows = emp::slice(TRANSITION_PROBS, ':');
        }

        if ((int)rows.size() != N_GENOTYPES) {
            std::cout << "Error: Not enough rows in transition matrix." << 
                " Attempting to assign " << rows.size() << 
                " to " << N_GENOTYPES << " genotypes. Rows: "
                << emp::to_string(rows) << std::endl;
            exit(1);
        }

        std::cout << "Transition probabilities: " << std::endl;
        for (int i = 0; i < N_GENOTYPES; i++) {
            emp::vector<std::string> sliced_probs = emp::slice(rows[i], ',');

            if ((int)sliced_probs.size() != N_GENOTYPES) {
                std::cout << "Error: Not enough columns in row." << 
                    " Attempting to assign " << sliced_probs.size() << 
                    " to " << N_GENOTYPES << " genotypes."
                    << std::endl;
                exit(1);
            }

            for (int j = 0; j < N_GENOTYPES; j++) {
                mut_rates[i][j] = emp::from_string<double>(sliced_probs[j]);
                std::cout << emp::from_string<double>(sliced_probs[j]) << " ";
            }
            std::cout << std::endl;
            
            if (mut_rates[i].GetWeight() != 1) {
                std::cout << "Error: Transition probabilities in row must sum to 1" << std::endl;
                exit(1);
            }
        }
        std::cout << std::endl;
    }

};
