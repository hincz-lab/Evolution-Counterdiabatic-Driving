# Evolution

Codes dealing with counter adiabatic control of evolutionary processes

## Running the model

First, clone this repository by executing the following command in a terminal:

```bash
git clone https://github.com/Peyara/Evolution.git
```

Currently, this repository contains a 1-dimensional
implementation of the model and the begginings of an
n-dimensional model (in progress).

The n-dimensional model has a dependency on the Empirical library. To run it, clone that as well:

```bash
git clone https://github.com/emilydolson/Empirical.git
git checkout memic_model # We're using a few features that haven't been merged into master yet
```

Then go into the directory for this repository:

```bash
cd Evolution
```

### 1-Dimensional Model

You can run the 1D model with the following commands:

```bash
make 1d  # compile the code
./1_dimension  # run the code
```

You will be prompted to enter the number of generations
to run each simulation for.

### N-dimensional Model

You can run the n-dimensional model with:

```bash
make nd  # compile the code
./n_dimensions  # run the code
```

This model has a few parameters:

- RANDOM_SEED (integer): the seed to initialize the random number generator with (for repeatability)
- GENERATIONS (integer): the number of generations to run for
- N_GENOTYPES (integer): the number of possible genotypes to mutate among (i.e. the number of dimensions to use for this realization of the model)
- K (floating point number): the carrying capacity
- DEATH_RATE (floating point number): The probability of an indiviudal dieing each generation (*d* in the equations)
- MAX_BIRTH_RATE (floating point number): The maximum possible birth rate (b0 in the equations)
- FITNESSES (string): Each genotype in the model needs a fitness. These fitnesses are expressed relative to a focal genotype (genotype 0), using the *s* parameter from the equations. Fitnesses can be specified directly on the command line or by providing the name of a file to look them up in. They should be specified as numbers separated by commas (see relative_fitnesses.dat for an example). The model needs a number of fitnesses equal to the number of genotypes. If it recieves one too few, it will assume the missing one belongs to the focal genotype and pre-pend a 0 to the list (since the focal genotype's fitness difference from itself is 0). Example: `-FITNESSES 0,.3,.2`.
- FITNESS_CHANGE_RULE (int): The goal of this model is to explore how we can steer evolution by changing relative fitnesses. There are a variety of different regimes that we might use: 0 (NONE) - relative fitnesses stay constant, 1 (VAR), 2 (VARCD) - a counteraadiabatic steering protocol. Options 1 and 2 currently only affect genotype 1's fitness, and ignore initial fitness values.
- INIT_POPS (string): Each genotype needs an initial population level. These are specified in the same way as FITNESSES (see init_pops.dat for an example).
- TRANSITION_PROBS (string): Each pair of genotypes needs probabilities for mutating between those genotypes (the probabilities can be different depending on the direction of the mutation). These probabilities are specified by providing a matrix. The diagonal of the matrix indicates the probabilities of not mutating. As with FITNESSES and INIT_POPS, this can either be provided in line or by referring to a file. The matrix is laid-out such that rows are the genotype being mutated from and columns are the genotype being mutated too. Within a row, numbers should be separated by commas. In a file rows can be separated with new-lines (example in transition_probs.dat). Because new-lines are hard to type as part of a command, rows should be separated with colons when the matrix is entered on the command line. Example: `-TRANSITION_PROBS .1,.9:.2,.8`. The matrix should be square (same number of rows and columns), and have a number of rows and columns equal to the number of genotypes. Numbers within a row should sum to 1 (since offspring need to have one of the available genotypes).

The values of these can be set with command line flags by placeing a dash before the name of the parameter you would like to modify and following it with the desired parameter value:

```bash
# For example, the following runs the model with a DEATH_RATE of .1 and GENERATIONS set to 100
./n_dimensions -DEATH_RATE .1 -GENERATIONS 100
```
