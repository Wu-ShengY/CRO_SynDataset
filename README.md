# Cyber Resilient OPF Dataset Synthesization

## Installation and usage

All models are implemented in Julia Language v.1.10.5 using [JuMP](https://github.com/jump-dev/JuMP.jl) modeling language for mathematical optimization and commercial [Mosek](https://github.com/MOSEK/Mosek.jl) or [Gurobi](https://github.com/jump-dev/Gurobi.jl) optimization solvers, which need to be licensed (free for academic use). 

The codes to implement CRO and CRO_Exp algorithms are placed in ```CRO_main.jl``` and ```CRO_EM_main.jl``` files, respectively. Make sure to active project environment using ```Project.toml``` and ```Manifest.toml``` located in the folder. 


### Virtual environment setup

First go into the CRO folder in the terminal

```cd your_file_path/CRO```

activate the julia environment by typing in ```julia```

Use the following command to setup the environments, then run the ```CRO_main.jl``` and ```CRO_EM_main.jl``` files with default arguments

```
]
activate .
instantiate
add PowerModels Statistics LinearAlgebra Distributions Random JuMP Gurobi DataFrames CSV Tables Plots LaTeXStrings CairoMakie StatsPlots ProgressMeter ArgParse Distributed Serialization
```

If there are some package conflicts, use ```Pkg.resolve()``` to resolve the conflicts

### REPL environment running

After ensuring that all packages are installed, type backspace to return to the julia environment

Run the scripts using ```include("CRO_main.jl")``` and ```include("CRO_EM_main.jl")```

### Terminal running

After environment setup, exit the REPL environment using

```exit()```

To run the CRO algorithm on the PJM 5-bus system, type the following command in the terminal:

```julia CRO_main.jl -a 50.0 -e 1.0```

which asks to compute the results for adjacency parameter 50.0MW and privacy loss 1.0. 

Similarly, to run the CRO_Exp algorithm on the IEEE 14-bus system, type:

```julia CRO_EM_main.jl -a 2.0 -e 1.0 -t 5```

which asks to compute the results from ```τ=1``` to ```τ=5``` , with adjacency parameter 2.0MW and privacy loss 1.0. 

---
