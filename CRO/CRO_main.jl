using PowerModels
using Statistics, LinearAlgebra, Distributions, Random
using JuMP, Gurobi, MosekTools, Mosek
using DataFrames, CSV, Tables, Plots
using LaTeXStrings
using CairoMakie, StatsPlots
using ProgressMeter, ArgParse
using Random
using Base.Threads
using Dates
using Distributed
include("aux_fun.jl")

# Load network data using PowerModels
cd(dirname(@__FILE__))
PowerModels.silence()

# net=load_network_data("./data/pglib_opf_case5_pjm.m")
net=load_network_data("./data/pglib_opf_case5_pjm__api.m")

# network parameters in compact format
net_c = Dict(:A => vcat(-ones(net[:N])' , net[:F] , -net[:F]), :B => vcat(ones(net[:N])' , -net[:F] , net[:F]), :c => vcat(ones(net[:N])'*net[:d] , -net[:F]*net[:d]-net[:f̅] , net[:F]*net[:d]-net[:f̅]))
# indexes for the nodes with loads
M_d= net[:d].!=0

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
            default = 50
        "--alpha", "-a"
            help = "Dataset Adjacency"
            arg_type = Float64
            default = 50
        "--loss", "-e"
            help = "privacy loss ε"
            arg_type = Float64
            default = 1.0
        "--magitude", "-m"
            help = "Magitudes of attack"
            arg_type = Float64
            default = 0.05
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
            default = 0.001
    end
    return parse_args(s)
end
args = parse_commandline()

att = Dict(:δ̅ => net[:d]*args["magitude"], :δ̲ => -net[:d]*args["magitude"], :Ψ =>args["Psi"], :τ=>args["tau"])
syn_set=Dict(:α =>args["alpha"] , :c_max =>0.0, :ϵ =>args["loss"], :γ=>args["gamma"], :β=>args["magitude"],  :s=>args["samples"])

# calculate the maximum cost of generators
syn_set[:c_max] = maximum(net[:c1][net[:gen_bus_ind]])

include("opt_fun.jl")


# To get numerical results of normal and post-attack OPF cost in PJM 5-bus systems
dict_pp, dict_cro=Results_CRO(net,net_c,syn_set[:s])
@show [mean(dict_pp[:C_OPF]) mean(dict_pp[:C_att_BO]) mean(dict_pp[:C_att_RO]) mean(dict_cro[:C_OPF]) mean(dict_cro[:C_att_BO])  mean(dict_cro[:C_att_RO]) ]

# Get results and serialize them
# CaseStudy_CRO(net,net_c,syn_set[:s])

# To get figures from saved results
# plot_CRO()

# To generate new figures from scratch
# plot_CRO_test(net,net_c,syn_set[:s])
