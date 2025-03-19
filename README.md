# Cyber Resilient OPF Dataset Synthesization

## Installation and usage

All models are implemented in Julia Language v.1.10.5 using [JuMP](https://github.com/jump-dev/JuMP.jl) modeling language for mathematical optimization and commercial [Mosek](https://github.com/MOSEK/Mosek.jl) and [Gurobi](https://github.com/jump-dev/Gurobi.jl) optimization solvers, which need to be licensed (free for academic use). 

The codes to implement CRO and CRO_Exp algorithms are placed in ```CRO_main.jl``` and ```CRO_EM_main.jl``` files, respectively. Make sure to active project environment using ```Project.toml``` and ```Manifest.toml``` located in the folder. 

First change the directory to the CRO folder:

```cd your_file_path/CRO```

### REPL environment

Use the following command to setup the environments, then run the ```CRO_main.jl``` and ```CRO_EM_main.jl``` files with default arguments

```
(CRO) pkg> activate .
(CRO) pkg> add Pkg
(CRO) pkg> instantiate .
```

### Terminal usage

To run the CRO algorithm on the PJM 5-bus system, type the following command in the terminal:

```julia CRO_main.jl -a 50.0 -e 1.0```

which asks to compute the results for adjacency parameter 50.0MW and privacy loss 1.0. 

Similarly, to run the CRO_Exp algorithm on the IEEE 14-bus system, type:

```julia CRO_EM_main.jl -a 2.0 -e 1.0 -t 5```

which asks to compute the results from ```τ=1``` to ```τ=5``` , with adjacency parameter 2.0MW and privacy loss 1.0. 

---
