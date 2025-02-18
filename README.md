# Cyber-Resilient-OPF-Dataset-Generation

## Installation and usage

All models are implemented in Julia Language v.1.8 using [JuMP](https://github.com/jump-dev/JuMP.jl) modeling language for mathematical optimization and commercial [Mosek](https://github.com/MOSEK/Mosek.jl) and [Gurobi](https://github.com/jump-dev/Gurobi.jl) optimization solvers, which need to be licensed (free for academic use). 

The codes to implement the two algorithms are placed in ```CRO_main.jl``` and ```CRO_EM_main.jl``` folders, respectively. Make sure to active project environment using ```Project.toml``` and ```Manifest.toml``` located in the folder. 

To run the WPO algorithm, ```cd``` to ```WPO``` and type the following command in the terminal:

```julia wpo_main.jl -a 15.0 -e 1.0```

which asks to compute the results for adjacency parameter 15.0 and privacy loss 1.0. 

Similarly, ```cd``` to ```TCO``` and type:

```julia tco_main.jl -a 15.0 -e 1.0```

to run the TCO algorithm. For more information on the settings, type

```julia tco_main.jl --help```

---
