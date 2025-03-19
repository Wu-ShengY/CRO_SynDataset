using PowerModels
using Statistics, LinearAlgebra, Distributions, Random
using JuMP, Gurobi
using DataFrames, CSV, Tables, Plots
using LaTeXStrings
using CairoMakie, StatsPlots
using ProgressMeter, ArgParse
using Distributed
using Serialization

include("aux_fun.jl")

# Load network data using PowerModels
cd(dirname(@__FILE__))
PowerModels.silence()

# net=load_network_data("./data/pglib_opf_case5_pjm.m")
# net=load_network_data("./data/pglib_opf_case5_pjm__api.m")
net=load_network_data("./data/pglib_opf_case14_ieee__api.m")

# network parameters in compact format
net_c = Dict(:A => vcat(-ones(net[:N])' , net[:F] , -net[:F]), :B => vcat(ones(net[:N])' , -net[:F] , net[:F]), :c => vcat(ones(net[:N])'*net[:d] , -net[:F]*net[:d]-net[:f̅] , net[:F]*net[:d]-net[:f̅]))
# indexes for the nodes with loads
M_d= net[:d].!=0

# Case Study Parameters
# 5bus: δ=0.05, Ψ=200, α=50
# 14bus: δ=0.15, Ψ=300, α=2
# 24bus: δ=0.15, Ψ=800, α=10
# 118bus: δ=0.15, Ψ=600, α=10

# parse arguments 
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--tau", "-t"
            help = "Number of Attacked Constraints"
            arg_type = Int64
            default = 5
        "--samples", "-s"
            help = "Number of Samples to run distribution"
            arg_type = Int64
            default = 10
        "--alpha", "-a"
            help = "Dataset Adjacency"
            arg_type = Float64
            default = 2
        "--loss", "-e"
            help = "privacy loss ε"
            arg_type = Float64
            default = 1.0
        "--magitude", "-m"
            help = "Magitudes of attack"
            arg_type = Float64
            default = 0.15
        "--Psi", "-p"
            help = "Parameters for constraint violation in OPF"
            arg_type = Float64
            default = 300
        "--beta", "-b"
            help = "Trade off parameters β in CRO algorithm"
            arg_type = Float64
            default = 5
        "--gamma", "-g"
            help = "regularization parameters γ in CRO algorithm"
            arg_type = Float64
            default = 0.0001
    end
    return parse_args(s)
end
args = parse_commandline()

att = Dict(:δ̅ => net[:d]*args["magitude"], :δ̲ => -net[:d]*args["magitude"], :Ψ =>args["Psi"], :τ=>args["tau"])
syn_set=Dict(:α =>args["alpha"] , :c_max =>0.0, :ϵ =>args["loss"], :γ=>args["gamma"], :β=>args["magitude"],  :s=>args["samples"])

# calculate the maximum cost of generators
syn_set[:c_max] = maximum(net[:c1][net[:gen_bus_ind]])

include("opt_fun.jl")


# To get attack damage on CRO_exp dataset under different τ
Results_CRO_EM(syn_set[:s])
