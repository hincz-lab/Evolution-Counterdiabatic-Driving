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
    VALUE(K, double, 10000, "Carrying capacity"),
    VALUE(DEATH_RATE, double, .05, "Death rate (d)"),
    VALUE(MAX_BIRTH_RATE, double, 2, "Maximum birth rate (b0)"),

    GROUP(PER_GENOTYPE_VALUES, "Per-genotype values"),
    VALUE(FITNESSES, std::string, "0,1", "Either a list of relative fitnesses, separated by commas, or a file containing them"),
    VALUE(FITNESS_CHANGE_RULE, int, 0, "Rule governing how fitnesses should change. 0 = NONE, 1 = VAR, 2 = VARCD"),
    VALUE(GENOTYPE_TO_SPEED_CONTROL, int, 0, "For fitness change rules that only apply to one genotype (currently all of them), which genotype should be changed?"),    
    VALUE(INIT_POPS, std::string, "100,10", "Either a list of initial population sizes, separated by commas, or a file containing them"),
    VALUE(TRANSITION_PROBS, std::string, ".95,.05:.05,.95", 
        "Either a matrix of transition probabilities or a file containing one. Rows are original genotype, columns are new one. Use commas to separate values within rows. In files, use newlines between rows. On command-line, use colons.")
)

enum class FITNESS_CHANGE_RULES { NONE=0, VAR=1, VARCD=2};

// Functions for changing s values over time
// Written by Shamreen

// Defining tanh based gen varying s function
double sVarCD(double x, double s)
{
//   double s;
  double ds,scd;
  s=(double)(0.00075+0.00075*tanh((x-500.)/270.));
  ds=0.00075/(270.*pow(cosh((x-500.)/270.),2.))/0.05;
  scd=s+(ds/pow((pow((0.0008-s),2.)+4.*0.0004*s),0.5));
  return scd;
}

// Defining tanh based gen varying s function
double sVar(double x, double s)
{
//   double s;
  s=(double)(0.00075+0.00075*tanh((x-500.)/270.));
  return s;
}

// Functions for extracting parameter lists from files

// Loads specified file into file object and prints appropriate warnings.
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

// Return the first line of the file. For use with one line files
std::string ExtractStringFromOneLineFile(std::string & fname, std::string param_name) {
    emp::File param_file = ExtractStringFromFile(fname, param_name);
    return param_file.front(); // Return first (theoretically only) line
}

// Returns the entire file. For use with multi-line files
emp::vector<std::string> ExtractMultiLineFile(std::string & fname, std::string param_name) {
    emp::File param_file = ExtractStringFromFile(fname, param_name);
    return param_file.GetAllLines(); // Return all lines
}


class NDimSim {
    private:
    emp::Random rnd; // Random-number generator
    emp::DataFile pop_sizes; // For tracking pops over time
    emp::DataFile pop_props; // For tracking genotype proportions of population over time
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
    int FITNESS_CHANGE_RULE;
    int GENOTYPE_TO_SPEED_CONTROL;
    std::string INIT_POPS;
    std::string TRANSITION_PROBS;

    public:
    NDimSim(EvoConfig & config) : pop_sizes("pop_sizes.csv"), pop_props("pop_props.csv") {
        Setup(config);
    }

    void Setup(EvoConfig & config) {
        rnd = emp::Random(config.RANDOM_SEED());

        // Localize config parameters
        // (we do this for efficiency)
        N_GENOTYPES = config.N_GENOTYPES();
        GENERATIONS = config.GENERATIONS();
        K = config.K();
        DEATH_RATE = config.DEATH_RATE();
        MAX_BIRTH_RATE = config.MAX_BIRTH_RATE();
        FITNESSES = config.FITNESSES();
        FITNESS_CHANGE_RULE = config.FITNESS_CHANGE_RULE();
        GENOTYPE_TO_SPEED_CONTROL = config.GENOTYPE_TO_SPEED_CONTROL();
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
        pop_props.AddVar(curr_gen,"generation");
        for (int i = 0; i < N_GENOTYPES; i++) {
            // Using a reference to vector contents is safe here
            // because vector size will never change
            pop_sizes.AddVar(current_pops[i],"pop" + emp::to_string(i));
            pop_props.AddFun((std::function<double()>)[i, this](){return current_pops[i]/emp::Sum(current_pops);},"prop" + emp::to_string(i));
        }
        pop_sizes.PrintHeaderKeys();
        pop_sizes.SetTimingRepeat(1);
        pop_props.PrintHeaderKeys();
        pop_props.SetTimingRepeat(1);


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

    // Updates fitnesses based on specified change rule
    void UpdateSs() {
        switch(FITNESS_CHANGE_RULE) {
            case (int)FITNESS_CHANGE_RULES::NONE:
                return;
                break;
            case (int)FITNESS_CHANGE_RULES::VAR:
                // This only really works for 1D right now
                rel_fitnesses[GENOTYPE_TO_SPEED_CONTROL] = sVar(curr_gen, rel_fitnesses[GENOTYPE_TO_SPEED_CONTROL]);
                break;
            case (int)FITNESS_CHANGE_RULES::VARCD:
                // This only really works for 1D right now
                rel_fitnesses[GENOTYPE_TO_SPEED_CONTROL] = sVarCD(curr_gen, rel_fitnesses[GENOTYPE_TO_SPEED_CONTROL]);
                break;
            default:
                std::cout << "Invalid fitness change rule. Defaulting to none." << std::endl;
                break;
        }
    }

    void RunStep() {
        UpdateSs(); // Update fitnesses as appropriate

        // Initialize new_pops to 0 so that we can accumulate the 
        // updated population sizes of each genotype
        for (size_t i = 0; i < new_pops.size(); i++) {
            new_pops[i] = 0;
        }

        // Handle birth/death/mutation
        for (int genotype = 0; genotype < N_GENOTYPES; genotype++) {
            // Store reference to current probability map, for simplicity
            emp::IndexMap & mut_probs = mut_rates[genotype];
            
            // Loop through all individuals of this genotype and figure out what happens to them
            for (int individual = 0; individual < current_pops[genotype]; individual++) {

                // Check if individual died
                if (rnd.GetDouble(1) - DEATH_RATE >= 0) {
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

    // Run for specified number of generations
    void Run() {
        for (int gen = 0; gen <= GENERATIONS; gen++) {
            curr_gen = gen;
            pop_sizes.Update(gen);
            pop_props.Update(gen);
            RunStep();
        }
    }

    // Pull fitness values out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeFitnesses() {
        if (emp::has_letter(FITNESSES)) {
            // FITNESSES is a file
            FITNESSES = ExtractStringFromOneLineFile(FITNESSES, "fitness");
        }

        // Split string up into smaller strings that contain individual values
        emp::vector<std::string> sliced_fits = emp::slice(FITNESSES, ',');

        // There should be N_GENOTYPES fitness values
        if ((int)sliced_fits.size() != N_GENOTYPES) {
            // If there are a different number of fitnesses, then
            // there's really nothing we can do
            std::cout << "Error: Not enough fitnesses supplied." << 
                " Attempting to assign " << sliced_fits.size() << 
                " to " << N_GENOTYPES << " genotypes." << std::endl;
            exit(1);
        }

        // Put the relative fitnesses we're using in the correct vector
        // and print them out
        std::cout << "Relative fitnesses: " << std::endl;
        for (int i = 0; i < N_GENOTYPES; i++) {
            rel_fitnesses[i] = emp::from_string<double>(sliced_fits[i]);
            std::cout << i << ": " << rel_fitnesses[i]  << " " << Birth(i) << std::endl;
        }
        std::cout << std::endl;
    }

    // Pull initial population values out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeInitPops() {
        if (emp::has_letter(INIT_POPS)) {
            // INIT_POPS is a filename
            INIT_POPS = ExtractStringFromOneLineFile(INIT_POPS, "initial population");
        }

        emp::vector<std::string> sliced_pops = emp::slice(INIT_POPS, ',');

        // Need a population size for each genotype
        if ((int)sliced_pops.size() != N_GENOTYPES) {
            std::cout << "Error: Not enough initial population sizes supplied." << 
                " Attempting to assign " << sliced_pops.size() << 
                " to " << N_GENOTYPES << " genotypes." << std::endl;
            exit(1);
        }

        // Put the initial population sizes we're using in the correct vector
        // and print them out
        std::cout << "Initial populations: " << std::endl;
        for (int i = 0; i < N_GENOTYPES; i++) {
            init_pops[i] = emp::from_string<double>(sliced_pops[i]);
            std::cout << i << ": " << init_pops[i] << std::endl;
        }
        std::cout << std::endl;
    }

    // Pull transition probabilities out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeTransitionProbs() {
        // First separate the rows
        emp::vector<std::string> rows;
        if (emp::has_one_of(TRANSITION_PROBS, "abcdefghijklmopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ")) {
            // INIT_POPS is a filename
            // This will technically fail if the file name is composed excelusively of 
            // lower case ns, numbers, and punctuation, but that seems unlikely.
            rows = ExtractMultiLineFile(TRANSITION_PROBS, "transition probabilities");
        } else {
            rows = emp::slice(TRANSITION_PROBS, ':');
        }

        // Need to have transition probabilities for everything
        if ((int)rows.size() != N_GENOTYPES) {
            std::cout << "Error: Not enough rows in transition matrix." << 
                " Attempting to assign " << rows.size() << 
                " to " << N_GENOTYPES << " genotypes. Rows: "
                << emp::to_string(rows) << std::endl;
            exit(1);
        }

        std::cout << "Transition probabilities: " << std::endl;
        for (int i = 0; i < N_GENOTYPES; i++) {
            // Now separate numbers within a row
            emp::vector<std::string> sliced_probs = emp::slice(rows[i], ',');

            // Need to have transition probabilities for everything
            if ((int)sliced_probs.size() != N_GENOTYPES) {
                std::cout << "Error: Not enough columns in row." << 
                    " Attempting to assign " << sliced_probs.size() << 
                    " to " << N_GENOTYPES << " genotypes."
                    << std::endl;
                exit(1);
            }

            // Fill in the mutation rates vector/IndexMap and print the
            // probabilities we're using
            for (int j = 0; j < N_GENOTYPES; j++) {
                mut_rates[i][j] = emp::from_string<double>(sliced_probs[j]);
                std::cout << emp::from_string<double>(sliced_probs[j]) << " ";
            }
            std::cout << std::endl;
            
            // Weights within each row have to sum to 1 because that's how
            // probability works
            if (abs(mut_rates[i].GetWeight() - 1) > .00000001) {
                std::cout << "Error: Transition probabilities in row must sum to 1. Sum is " << mut_rates[i].GetWeight() << std::endl;
                exit(1);
            }
        }
        std::cout << std::endl;
    }

};
