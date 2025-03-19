include("aux_fun.jl")
include("opt_fun.jl")

# Load network data using PowerModels
cd(dirname(@__FILE__))
PowerModels.silence()

# net=load_network_data("./data/pglib_opf_case5_pjm.m")
net=load_network_data("./data/pglib_opf_case5_pjm__api.m")

# parse arguments 
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--tau", "-t"
            help = "Number of Attacked Constraints"
            arg_type = Int64
            default = 10
        "--samples", "-s"
            help = "Number of Samples to run distribution"
            arg_type = Int64
            default = 50
        "--alpha", "-a"
            help = "Dataset Adjacency"
            arg_type = Float64
            default = 5.0
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
            default = 200
        "--beta", "-b"
            help = "Trade off parameters β in CRO algorithm"
            arg_type = Float64
            default = 5
        "--gamma", "-g"
            help = "regularization parameters γ in CRO algorithm"
            arg_type = Float64
            default = 0.00001
    end
    return parse_args(s)
end
args = parse_commandline()

att = Dict(:δ̅ => net[:d]*args["magitude"], :δ̲ => -net[:d]*args["magitude"], :Ψ =>args["Psi"], :τ=>args["tau"])
syn_set=Dict(:α =>args["alpha"] , :c_max =>0.0, :ϵ =>args["loss"], :γ=>args["gamma"], :β=>args["magitude"],  :s=>args["samples"])

# To get numerical results of normal and post-attack OPF cost in PJM 5-bus systems
# dict_pp, dict_cro=Results_CRO(100)
# @show [mean(dict_cro[:C̃]) std(dict_cro[:C̃]) mean(dict_cro[:C_OPF]) std(dict_cro[:C_OPF]) mean(dict_cro[:C_att_BO]) std(dict_cro[:C_att_BO]) mean(dict_cro[:C_att_RO]) std(dict_cro[:C_att_RO])]

# Get results and serialize them
# CaseStudy_CRO(syn_set[:s])

# To get figures from saved results
# plot_CRO()

# To generate new figures from scratch
plot_CRO_test(syn_set[:s])
