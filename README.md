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

The values of these can be set with command line flags by placeing a dash before the name of the parameter you would like to modify and following it with the desired parameter value:

```bash
# For example, the following runs the model with a DEATH_RATE of .1 and GENERATIONS set to 100
./n_dimensions -DEATH_RATE .1 -GENERATIONS 100
```
