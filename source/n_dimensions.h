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

    GROUP(FITNESS_CHANGE_PARAMETERS, "Parameters associated with various fitness change rules"),
    VALUE(FITNESS_CHANGE_RULE, int, 0, "Rule governing how fitnesses should change. 0 = NONE, 1 = VAR, 2 = VARCD, 3 = Drug with increasing dose, 4 = Drug with fixed dose, 5 = CD Driving prescription specified in file."),
    VALUE(GENOTYPE_TO_DRIVE, int, 0, "For fitness change rules that only apply to one genotype (VAR and VARCD), which genotype should be changed?"),    
    VALUE(TIME_STEPS_BEFORE_RAMP_UP, int, 0, "For fitness change rule 3, how long to wait before we start increasing concentration"),
    VALUE(DRUG_DOSE, double, .00015, "For fitness change rules 3 and 4"),
    VALUE(CD_DRIVING_PRESCRIPTION, std::string, "driving.csv", "File containing driving prescription for use with fitness change rule 5"),

    GROUP(PER_GENOTYPE_VALUES, "Per-genotype values"),
    VALUE(FITNESSES, std::string, "0,1", "Either a list of relative fitnesses, separated by commas, or a file containing them. These are the starting ftnesses."),
    VALUE(IC50S, std::string, "-6.0,-5.0", "For environments simulating the application of a drug, what are the IC50 values for each genotype? Specify as list of values or name of file containing them."),
    VALUE(G_DRUGLESSES, std::string, "1,1", "For environments simulating the application of a drug, what are the growth rates in absence of the drug? Specify as list of values or name of file containing them."),
    VALUE(CS, std::string, "1,1", "Constants describing the shape of the hill functions relating dose to fitness for each genotype."),
    VALUE(INIT_POPS, std::string, "100,10", "Either a list of initial population sizes, separated by commas, or a file containing them"),
    VALUE(TRANSITION_PROBS, std::string, ".95,.05:.05,.95", 
        "Either a matrix of transition probabilities or a file containing one. Rows are original genotype, columns are new one. Use commas to separate values within rows. In files, use newlines between rows. On command-line, use colons.")
)

enum class FITNESS_CHANGE_RULES { NONE=0, VAR=1, VARCD=2, INCREASING_DRUG=3, CONSTANT_DRUG=4, CD_PRESCRIPTION=5};

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
    emp::vector<long double> rel_fitnesses;
    emp::vector<long double> init_pops;
    emp::vector<long double> current_pops;
    emp::vector<long double> new_pops;
    emp::vector<emp::IndexMap> mut_rates;

    // Drug environment specific per-genotypes values
    emp::vector<long double> drugless_fitnesses;
    emp::vector<long double> cs;
    emp::vector<long double> IC50s;
    long double max_fit;

    // CD driving prescription (if necessary)
    emp::vector<emp::vector<long double>> cd_prescription_data;

    // Localized config parameters
    int N_GENOTYPES;
    int GENERATIONS;
    double K;
    double DEATH_RATE;
    double MAX_BIRTH_RATE;
    std::string FITNESSES;
    int FITNESS_CHANGE_RULE;
    int TIME_STEPS_BEFORE_RAMP_UP;
    double DRUG_DOSE;
    int GENOTYPE_TO_DRIVE;
    std::string INIT_POPS;
    std::string TRANSITION_PROBS;
    std::string IC50S;
    std::string CS;
    std::string G_DRUGLESSES;
    std::string CD_DRIVING_PRESCRIPTION;

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
        DRUG_DOSE = config.DRUG_DOSE();
        GENOTYPE_TO_DRIVE = config.GENOTYPE_TO_DRIVE();
        INIT_POPS = config.INIT_POPS();
        TRANSITION_PROBS = config.TRANSITION_PROBS();
        IC50S = config.IC50S();
        CS = config.CS();        
        G_DRUGLESSES = config.G_DRUGLESSES();        
        TIME_STEPS_BEFORE_RAMP_UP = config.TIME_STEPS_BEFORE_RAMP_UP();        
        CD_DRIVING_PRESCRIPTION = config.CD_DRIVING_PRESCRIPTION();        

        // Set-up per-genotype values

        // Fitnesses (s in the equation) are relative to first genotype 
        // The first genotype's s will be 0
        InitializeFitnesses();

        // Initial population sizes of each genotype
        InitializeInitPops();

        // Initialize IC50 values for each genotype
        InitializeIC50s();

        // Initialize frugless growth rates for each genotype
        InitializeGDruglesses();

        // Initialize c values for each genotype
        InitializeCs();

        if (FITNESS_CHANGE_RULE == (int)FITNESS_CHANGE_RULES::CONSTANT_DRUG) {
            sDrugConcentration(DRUG_DOSE);
        } else if (FITNESS_CHANGE_RULE == (int)FITNESS_CHANGE_RULES::CD_PRESCRIPTION) { 
            cd_prescription_data = emp::File(CD_DRIVING_PRESCRIPTION).ToData<long double>();
        }

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

    /* Birth rate of other genotypes (gets called for 
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
            case (int)FITNESS_CHANGE_RULES::CONSTANT_DRUG:
                return;
                break;
            case (int)FITNESS_CHANGE_RULES::VAR:
                // This only really works for 1D
                rel_fitnesses[GENOTYPE_TO_DRIVE] = sVar(curr_gen, rel_fitnesses[GENOTYPE_TO_DRIVE]);
                break;
            case (int)FITNESS_CHANGE_RULES::VARCD:
                // This only really works for 1D
                rel_fitnesses[GENOTYPE_TO_DRIVE] = sVarCD(curr_gen, rel_fitnesses[GENOTYPE_TO_DRIVE]);
                break;
            case (int)FITNESS_CHANGE_RULES::INCREASING_DRUG:
                sDrugIncrease(curr_gen);
                break;
            case (int)FITNESS_CHANGE_RULES::CD_PRESCRIPTION:
                rel_fitnesses = cd_prescription_data[curr_gen];
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


    // Methods for changing s values over time

    // Set s for a sepcific drug concentration
    void sDrugConcentration(double concentration) {
        if (concentration > 0) {
            for (int genotype = 0; genotype < N_GENOTYPES; genotype++) {
                rel_fitnesses[genotype] = drugless_fitnesses[genotype]/(1 + exp((IC50s[genotype] - emp::Log10(concentration))/cs[genotype]));
            }
        } else {
            for (int genotype = 0; genotype < N_GENOTYPES; genotype++) {
                rel_fitnesses[genotype] = drugless_fitnesses[genotype];
            }
        }
        // std::cout << emp::to_string(rel_fitnesses) << std::endl;
        for (int genotype = 0; genotype < N_GENOTYPES; genotype++) {
            rel_fitnesses[genotype] = rel_fitnesses[15]/rel_fitnesses[genotype] - 1;
        }
    }

    void sDrugIncrease(double t) {
        long double concentration = 0;
        if (t >= TIME_STEPS_BEFORE_RAMP_UP/2) {
            concentration = DRUG_DOSE/(1 + exp(-.002*(t - TIME_STEPS_BEFORE_RAMP_UP - 110)));
        }
        sDrugConcentration(concentration);
   }

    // void sDrugIncrease(double t) {
    //     // This function ramps up over 2M generations
    //     // double concentration = tanh(t/500000)*10000;
        
    //     // This one ramps up over 100 generations
    //     // double concentration = tanh(.03*t);
    //     // long double concentration = .1/(1 + exp(-.224*(t - 75)));
    //     long double concentration = 0;
    //     if (t >= TIME_STEPS_BEFORE_RAMP_UP) {
    //         concentration = .01/(1 + exp(-.07*(t - TIME_STEPS_BEFORE_RAMP_UP - 110)));
    //     }
    //     // std::cout << "Concentration: " << concentration << std::endl;
    //     sDrugConcentration(concentration);
    // }

    // The following two functions were written by Shamreen Iram
    // Defining tanh based gen varying s function
    double sVarCD(double x, double s) {
        //   double s;
        double ds,scd;
        s=(double)(0.00075+0.00075*tanh((x-500.)/270.));
        ds=0.00075/(270.*pow(cosh((x-500.)/270.),2.))/0.05;
        scd=s+(ds/pow((pow((0.0008-s),2.)+4.*0.0004*s),0.5));
        return scd;
    }

    // Defining tanh based gen varying s function
    double sVar(double x, double s) {
        //   double s;
        s=(double)(0.00075+0.00075*tanh((x-500.)/270.));
        return s;
    }

    emp::vector<emp::IndexMap> GetMutRates() {return mut_rates;}

    emp::vector<long double> ExtractVectorFromConfig(std::string param, std::string name, std::string plural) {
        emp::vector<long double> result(N_GENOTYPES);

        if (emp::has_letter(param)) {
            // param is a file
            param = ExtractStringFromOneLineFile(param, name);
        }

        // Split string up into smaller strings that contain individual values
        emp::vector<std::string> sliced_param = emp::slice(param, ',');

        // There should be N_GENOTYPES values
        if ((int)sliced_param.size() != N_GENOTYPES) {
            // If there are a different number of then
            // there's really nothing we can do
            std::cout << "Error: Not enough " << plural << "supplied." << 
                " Attempting to assign " << sliced_param.size() << 
                " to " << N_GENOTYPES << " genotypes." << std::endl;
            exit(1);
        }

        // Put the values we're using in the correct vector
        // and print them out
        std::cout << plural << ": " << std::endl;
        for (int i = 0; i < N_GENOTYPES; i++) {
            result[i] = emp::from_string<long double>(sliced_param[i]);
            std::cout << i << ": " << result[i]  << " " << std::endl;
        }
        std::cout << std::endl;
        return result;
    }

    // Pull fitness values out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeFitnesses() {
        rel_fitnesses = ExtractVectorFromConfig(FITNESSES, "fitness", "fitnesses");
    }

    // Pull initial population values out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeInitPops() {
        init_pops = ExtractVectorFromConfig(INIT_POPS, "initial population", "initial populations");
    }

    // Pull IC50 values out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeIC50s() {
        IC50s = ExtractVectorFromConfig(IC50S, "IC50", "IC50s");
    }

    // Pull drugless growth rates out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeGDruglesses() {
        drugless_fitnesses = ExtractVectorFromConfig(G_DRUGLESSES, "g_drugless", "g_druglesses");
        max_fit = emp::FindMax(drugless_fitnesses);
    }

    // Pull c values out of config parameter and put them in the
    // appropriate vector, giveing warnings and errors as neccessary.
    void InitializeCs() {
        cs = ExtractVectorFromConfig(CS, "C", "Cs");
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
