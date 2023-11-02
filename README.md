# sEMTor.jl

Julia simulation for epithelial-to-mesenchymal transitions (EMT).

## Stand-alone EMT simulator (sEMTor)

For trying the model, one can use the JavaScript version of the simulator:
[semtor.github.io](https://semtor.github.io).

The JavaScript variant implements the same model, but without abilities for data collection. 
For running scientific simulations and performing parameter studies, please use the Julia version.

## Installation

The simulator is written in Julia, and can be installed as a Julia package.
To install Julia, please follow the instructions on the [Julia website](https://julialang.org/downloads/).

To install the simulator, start Julia and run the following commands:

```julia
using Pkg
Pkg.add("github.com/SteffenPL/sEMTor.jl")  # or add 'github.com/SteffenPL/sEMTor.jl' in pkg mode
``` 

## Usage

The parameters for the EMT model can be provided either as TOML files or as XLSX files. 
Here, we show how to load XLSX files provided in the repository folder.

```julia
using SEMTor 

# Load parameters from XLSX files into a dictionary
params = load_parameters("params.xlsx")

# Modify parameters (if needed)
setvalue!(params, "emt.lifespan.min", 9.0)
setvalue!(params, "emt.lifespan.max", 10.0)

# Convert to a static type for efficient simulations
p = StaticParameters(params)

states = simulate(p)

# Plot the results
plot_states(states)
```

### Statistical analysis

The simulator can be used to perform statistical analysis of the EMT model. 
For this is it useful to run simulations in parallel. This requires to start Julia with multiple threads:

```bash
julia -t auto
```

The following example shows how to run 100 simulations in parallel and collect the results in a `DataFrame`:

```julia
using DataFrames
using SEMTOR 

# Load parameters from XLSX files into a dictionary
p = StaticParameters(load_parameters("params.xlsx"))

# Run 100 simulations in parallel
ensemble_sims = simulate_parallel(p, 100)

# Calculate statistics from ensemble simulations
statistics = compute_statistics(ensemble_sims, p)

# Save different types of statistics in separate files
save_statistics("output_folder", statistics)
```

### Heterogeneous input parameters

In the XLSX file and the TOML file, it is possible to define parameter ranges instead of single values.
This allows to run simulations with heterogeneous input parameters.

Ranges can be defined as follows, where 'a', 'b' and 'c' need to be replaced with numbers and 'p' a percentage value:
```
(a, b)  # Uniform distribution between a and b
(a, b, p% c)  # Uniform distribution between a and b, with p% of the values are set to 'c'
```

We note that values can be set to 'Inf' for infinity.