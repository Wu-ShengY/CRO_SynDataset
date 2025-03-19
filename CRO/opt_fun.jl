# standard DC-OPF
function OPF(net, d)
    # original DC-OPF problem
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 0)
    # model variables
    @variable(model, p[1:net[:N]])
    # model objective
    @objective(model, Min, net[:c1]'p)
    # OPF equations
    @constraint(model, λ, ones(net[:N])' * (p .- d) .>= 0)
    @constraint(model, μ̅, net[:f̅] .>= net[:F] * (p .- d))
    @constraint(model, μ̲, net[:f̅] .>= -net[:F] * (p .- d))
    @constraint(model, γ̅, -p .>= -net[:p̅])
    @constraint(model, γ̲, -net[:p̲] .>= -p)
    # solve model
    optimize!(model)
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        sol = Dict(:status => "$(termination_status(model))",
            :obj => JuMP.objective_value(model),
            :p => JuMP.value.(p),
            :CPUtime => solve_time(model)
        )
    else
        sol = Dict(:status => termination_status(model))
    end
    return sol
end

# standard DC-OPF in compact form (min c'p  s.t. Ap+Bd+C ⩽ 0)
function DC_OPF_compact_notation(net_c)
    K = size(net_c[:A], 1)
    # original DC-OPF problem
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_silent(model)
    set_optimizer_attribute(model, "OutputFlag", 0)
    # model variables
    @variable(model, p[1:net[:N]])
    # model objective
    @objective(model, Min, net[:c1]'p)
    # OPF equations
    @constraint(model, con[k=1:K], net_c[:A][k, :]' * p + net_c[:c][k] .<= 0)
    @constraint(model, γ̅, -p .>= -net[:p̅])
    @constraint(model, γ̲, -net[:p̲] .>= -p)
    # solve model
    optimize!(model)
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        sol = Dict(:status => "$(termination_status(model))",
            :obj => JuMP.objective_value(model),
            :p => JuMP.value.(p),
            :CPUtime => solve_time(model)
        )
    else
        sol = Dict(:status => termination_status(model))
    end
    return sol
end

# compact DC-OPF with slack variable on the line flow constraints to ensure feasibility
function DC_OPF_compact_feas(net_c)
    K = size(net_c[:A], 1)
    # original DC-OPF problem
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 0)
    # generation variables
    @variable(model, p[1:net[:N]])
    # slack variables
    @variable(model, v[1:net[:E]] .>= 0)
    @variable(model, m[1:K])
    @constraint(model, m[1] == 0)
    @constraint(model, m[2:1+net[:E]] .== v)
    @constraint(model, m[2+net[:E]:K] .== v)
    # model objective
    @objective(model, Min, net[:c1]'p + att[:Ψ] * (ones(net[:E])' * v))
    # OPF equations
    @constraint(model, con[k=1:K], net_c[:A][k, :]' * p + net_c[:c][k] .- m[k] .<= 0)
    @constraint(model, γ̅, -p .>= -net[:p̅])
    @constraint(model, γ̲, -net[:p̲] .>= -p)
    # solve model
    optimize!(model)
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        sol = Dict(:status => "$(termination_status(model))",
            :obj => JuMP.objective_value(model),
            :p => JuMP.value.(p),
            :CPUtime => solve_time(model),
            :v => JuMP.value.(v)
        )
    else
        sol = Dict(:status => termination_status(model))
    end
    return sol
end

# load redistribution attack model in bilevel optimization settings 
function AttackModel_BO(net, att)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 0)
    # upper-level variables 
    @variable(model, δ[1:net[:N]])
    # primal lower-level variables
    @variable(model, p[1:net[:N]])
    # dual variables
    @variable(model, λ)
    @variable(model, μ̅[1:net[:E]] >= 0)
    @variable(model, μ̲[1:net[:E]] >= 0)
    @variable(model, γ̅[1:net[:N]] >= 0)
    @variable(model, γ̲[1:net[:N]] >= 0)
    # upper-level objective function 
    @objective(model, Max, net[:c1]'p)
    # attack vector constraints 
    @constraint(model, δ .<= att[:δ̅])
    @constraint(model, δ .>= att[:δ̲])
    @constraint(model, ones(net[:N])'δ == 0)
    # KKT conditions of OPF
    # primal constriants 
    @constraint(model, ones(net[:N])' * (p .- net[:d] .- δ) .== 0)
    @constraint(model, -net[:F] * (p .- net[:d] .- δ) .>= -net[:f̅])
    @constraint(model, net[:f̅] .>= -net[:F] * (p .- net[:d] .- δ))
    @constraint(model, -p .>= -net[:p̅])
    @constraint(model, -net[:p̲] .>= -p)

    # stationary constraints
    @constraint(model, net[:c1] + λ * ones(net[:N]) + net[:F]' * μ̅ - net[:F]' * μ̲ + γ̅ - γ̲ .== 0)
    # # strong duality 
    # complementarity conditions formulated as SOS1 constraints
    @variable(model, aux_μ̅[1:net[:E]] .>= 0)
    @variable(model, aux_μ̲[1:net[:E]] .>= 0)
    @variable(model, aux_γ̅[1:net[:E]] .>= 0)
    @variable(model, aux_γ̲[1:net[:E]] .>= 0)
    for i in 1:net[:E]
        @constraint(model, aux_μ̅[i] .== net[:f̅][i] .- net[:F][i, :]' * (p .- net[:d] .- δ))
        @constraint(model, [μ̅[i], aux_μ̅[i]] in SOS1())

        @constraint(model, aux_μ̲[i] .== net[:F][i, :]' * (p .- net[:d] .- δ) .+ net[:f̅][i])
        @constraint(model, [μ̲[i], aux_μ̲[i]] in SOS1())
    end
    for i in 1:net[:N]
        @constraint(model, aux_γ̅[i] .== net[:p̅][i] .- p[i])
        @constraint(model, [γ̅[i], aux_γ̅[i]] in SOS1())

        @constraint(model, aux_γ̲[i] .== p[i] .- net[:p̲][i])
        @constraint(model, [γ̲[i], aux_γ̲[i]] in SOS1())
    end
    # solve primal-dual model
    optimize!(model)
    @info("done solving the Attack_KKT problem: $(termination_status(model))")
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        sol = Dict(:status => termination_status(model),
            :obj => JuMP.objective_value(model),
            :p => JuMP.value.(p),
            :CPUtime => solve_time(model),
            :δ => JuMP.value.(δ),
        )
    else
        sol = Dict(:status => termination_status(model))
    end
end

# load redistribution attack model in bilevel optimization settings allowing constraint violations
function AttackModel_BO_feas(net,att)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 0)
    # upper-level variables 
    @variable(model, δ[1:net[:N]])
    # primal lower-level variables
    @variable(model, p[1:net[:N]])
    @variable(model, v[1:net[:E]].>=0)
    # dual variables
    @variable(model, λ)
    @variable(model, θ[1:net[:E]]>=0) # dual for v>=0
    @variable(model, μ̅[1:net[:E]]>=0)
    @variable(model, μ̲[1:net[:E]]>=0)
    @variable(model, γ̅[1:net[:N]]>=0)
    @variable(model, γ̲[1:net[:N]]>=0)
    # upper-level objective function 
    @objective(model, Max, net[:c1]'p+ att[:Ψ]*(ones(net[:E])'*v))
    # attack vector constraints 
    @constraint(model, δ.<=att[:δ̅])
    @constraint(model, δ.>=att[:δ̲])
    @constraint(model, ones(net[:N])'δ==0)
    # KKT conditions of OPF
    # primal constriants 
    @constraint(model, ones(net[:N])'*(p .- net[:d].-δ) .== 0)
    @constraint(model, - net[:F]*(p .- net[:d].-δ) .>= - net[:f̅] .-v)
    @constraint(model, net[:f̅].+v .>= - net[:F]*(p .- net[:d].-δ))
    @constraint(model, -p .>= -net[:p̅])
    @constraint(model, -net[:p̲] .>= -p)

    # stationary constraints
    @constraint(model, net[:c1] + λ * ones(net[:N]) + net[:F]' * μ̅ - net[:F]' * μ̲ + γ̅ - γ̲ .== 0)
    @constraint(model, att[:Ψ]*(ones(net[:E])) .- μ̅ .- μ̲ .- θ== 0)
    # # strong duality 
    # complementarity conditions formulated as SOS1 constraints
    @variable(model,aux_μ̅[1:net[:E]].>=0)
    @variable(model,aux_μ̲[1:net[:E]].>=0)
    @variable(model,aux_γ̅[1:net[:E]].>=0)
    @variable(model,aux_γ̲[1:net[:E]].>=0)
    @variable(model,aux_θ[1:net[:E]].>=0)
    for i in 1:net[:E]
        @constraint(model, aux_μ̅[i] .== net[:f̅][i].+v[i] .- net[:F][i,:]'*(p .- net[:d].-δ))
        @constraint(model, [μ̅[i] , aux_μ̅[i]] in SOS1())

        @constraint(model, aux_μ̲[i] .== net[:F][i,:]'*(p .- net[:d].-δ) .+ net[:f̅][i].+v[i])
        @constraint(model, [μ̲[i] , aux_μ̲[i]] in SOS1())
        @constraint(model, aux_θ[i] .== v[i])
        @constraint(model, [θ[i] , aux_θ[i]] in SOS1())
    end
    for i in 1:net[:N]
        @constraint(model, aux_γ̅[i] .== net[:p̅][i] .- p[i])
        @constraint(model, [γ̅[i], aux_γ̅[i]] in SOS1())

        @constraint(model, aux_γ̲[i] .== p[i] .- net[:p̲][i])
        @constraint(model, [γ̲[i] , aux_γ̲[i]] in SOS1())
    end
    # solve primal-dual model
    optimize!(model)
    @info("done solving the Attack_BO_feas problem: $(termination_status(model))")
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        sol = Dict(:status => termination_status(model),
        :obj => JuMP.objective_value(model),
        :p => JuMP.value.(p),
        :CPUtime => solve_time(model),
        :δ => JuMP.value.(δ),
        )
    else
        sol = Dict(:status => termination_status(model))
    end
end

# load redistribution attack model in robust optimization settings 
function AttackModel_RO(net, net_c, att)
    K = size(net_c[:A], 1)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 0)

    # primal lower-level variables
    @variable(model, p[1:net[:N]])
    # slack variable for line capacity
    @variable(model, v[1:net[:E]] .>= 0)
    # OPF objective function 
    @objective(model, Min, net[:c1]'p + att[:Ψ] * (ones(net[:E])' * v))
    # dual variables 
    @variable(model, μ̅[1:K, 1:net[:N]])
    @variable(model, μ̲[1:K, 1:net[:N]])
    @variable(model, λ[1:K])
    # primal problem constraint 
    @constraint(model, -net_c[:A][1, :]'p .- net_c[:c][1, :] .- μ̅[1, :]' * att[:δ̅] .+ μ̲[1, :]' * att[:δ̲] .>= 0)
    @constraint(model, θk1[k=2:net[:E]+1], -net_c[:A][k, :]'p .- net_c[:c][k, :] .- μ̅[k, :]' * att[:δ̅] .+ μ̲[k, :]' * att[:δ̲] .+ v[k-1] .>= 0)
    @constraint(model, θk2[k=net[:E]+2:K], -net_c[:A][k, :]'p .- net_c[:c][k, :] .- μ̅[k, :]' * att[:δ̅] .+ μ̲[k, :]' * att[:δ̲] .+ v[k-net[:E]-1] .>= 0)
    # stationarity constraint 
    @constraint(model, ζ[k=1:K], net_c[:B][k, :] .- μ̅[k, :] .+ μ̲[k, :] .- ones(net[:N]) .* λ[k] .== 0)
    # power generation constraints
    @constraint(model, ω̅, -p .>= -net[:p̅])
    @constraint(model, ω̲, p .>= net[:p̲])
    # primal feasibility
    @constraint(model, π̅[k=1:K, n=1:net[:N]], μ̅[k, n] .>= 0)
    @constraint(model, π̲[k=1:K, n=1:net[:N]], μ̲[k, n] .>= 0)
    # solve primal-dual model
    optimize!(model)
    @info("done solving the Feasible AttackModel problem: $(termination_status(model))")
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        # get the value of dual variables 
        ζ_matrix = zeros(K, net[:N])
        π̅_matrix = zeros(K, net[:N])
        π̲_matrix = zeros(K, net[:N])
        v_matrix = zeros(net[:E])
        for k in 1:K
            ζ_matrix[k, :] = JuMP.dual.(ζ[k])
            for n in 1:net[:N]
                π̅_matrix[k, n] = JuMP.dual.(π̅[k, n])
                π̲_matrix[k, n] = JuMP.dual.(π̲[k, n])
            end
            for l in 1:net[:E]
                v_matrix[l] = JuMP.value(v[l])
            end
        end
        #calculate the attack vector
        sol = Dict(:status => status,
            :obj => JuMP.objective_value(model),
            :p => JuMP.value.(p),
            :CPUtime => solve_time(model),
            :ζ => ζ_matrix,
            :ω̅ => JuMP.dual.(ω̅),
            :ω̲ => JuMP.dual.(ω̲),
            :π̅ => π̅_matrix,
            :π̲ => π̲_matrix,
            :v => v_matrix
        )
    else
        sol = Dict(:status => status)
    end
end

# RO attacks that only attack constraints in the set W_Set 
function AttackModel_RO_k(net,net_c,att,W_Set)
    K = size(net_c[:A],1)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer, "OutputFlag"=>0))

    # primal lower-level variables
    @variable(model, p[1:net[:N]])
    # slack variable for line capacity
    @variable(model, v[1:net[:E]].>=0)
    # OPF objective function 
    @objective(model, Min, net[:c1]'p+ att[:Ψ]*(ones(net[:E])'*v))
    # auxiliary variable for compact power balance and line capacity constraint
    @variable(model, m[1:K].>=0)
    @constraint(model, m[1].==0)
    @constraint(model, m[2:net[:E]+1].==v)
    @constraint(model, m[net[:E]+2:K].==v)
    # dual variables 
    @variable(model, μ̅[1:K,1:net[:N]].>=0)
    @variable(model, μ̲[1:K,1:net[:N]].>=0)
    @variable(model, λ[1:K])
    for k in 1:K
        if k in W_Set
            # primal problem constraint
            @constraint(model, -net_c[:A][k,:]'p .- net_c[:c][k,:] .- μ̅[k,:]'*att[:δ̅] .+ μ̲[k,:]'*att[:δ̲] .+ m[k] .>= 0)
             # stationarity constraint 
            @constraint(model, net_c[:B][k,:] .- μ̅[k,:] .+ μ̲[k,:] .- ones(net[:N]) .* λ[k] .== 0)
        else
            @constraint(model, -net_c[:A][k,:]'p .- net_c[:c][k,:] .+ m[k] .>= 0)
        end
    end
    # power generation constraints
    @constraint(model, ω̅,-p .>= -net[:p̅])
    @constraint(model, ω̲,p .>= net[:p̲])

    # solve primal-dual model
    optimize!(model)    
    # @info("done solving the Feasible Attack_k problem: $(termination_status(model))")
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        v_matrix = zeros(net[:E])
        for l in 1:net[:E]
            v_matrix[l]=JuMP.value(v[l])
        end
        #calculate the attack vector
        sol = Dict(:status => status,
        :obj => JuMP.objective_value(model),
        :p => JuMP.value.(p),
        :CPUtime => solve_time(model),
        :v => JuMP.value.(v_matrix)
        )
    else
        sol = Dict(:status => status)
    end
    return sol
end

# standard post processing algorithm on the obsfuscated data to generate synthetic load datasets
function OPF_PP_Demand(net, d̃_0, C̃_d)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 0)
    # upper-level variables 
    @variable(model, d̃[1:net[:N]] .>= 0)
    @variable(model, r_d)
    @variable(model, r_c)
    # primal lower-level variables
    @variable(model, p[1:net[:N]])
    # dual variables
    @variable(model, λ >= 0)
    @variable(model, μ̅[1:net[:E]] .>= 0)
    @variable(model, μ̲[1:net[:E]] .>= 0)
    @variable(model, γ̅[1:net[:N]] .>= 0)
    @variable(model, γ̲[1:net[:N]] .>= 0)
    # objective function
    @objective(model, Min, syn_set[:γ] * r_d + r_c)
    # objective aux constraints 
    @constraint(model, [r_d; M_d .* d̃ .- d̃_0] in MOI.NormOneCone(1 + length(d̃)))
    @constraint(model, aux_con, [r_c; C̃_d .- net[:c1]'p] in MOI.NormOneCone(2))
    # primal constriants 
    @constraint(model, ones(net[:N])' * (p .- M_d .* d̃) .>= 0)
    @constraint(model, -net[:F] * (p .- M_d .* d̃) .>= -net[:f̅])
    @constraint(model, net[:f̅] .>= -net[:F] * (p .- M_d .* d̃))
    @constraint(model, p .<= net[:p̅])
    @constraint(model, p .>= net[:p̲])
    # dual constraints
    @constraint(model, λ * ones(net[:N]) - net[:F]' * μ̅ + net[:F]' * μ̲ - γ̅ + γ̲ .== net[:c1])
    # # strong duality 
    # complementarity conditions formulated as SOS1 constraints
    @variable(model, aux_μ̅[1:net[:E]] .>= 0)
    @variable(model, aux_μ̲[1:net[:E]] .>= 0)
    @variable(model, aux_γ̅[1:net[:N]] .>= 0)
    @variable(model, aux_γ̲[1:net[:N]] .>= 0)
    @variable(model, aux_λ >= 0)

    @constraint(model, aux_λ == ones(net[:N])' * (p .- M_d .* d̃))
    @constraint(model, [λ, aux_λ] in SOS1())
    for i in 1:net[:E]
        @constraint(model, aux_μ̅[i] == net[:f̅][i] - net[:F][i, :]' * (p .- M_d .* d̃))
        @constraint(model, [μ̅[i], aux_μ̅[i]] in SOS1())
        @constraint(model, aux_μ̲[i] == net[:F][i, :]' * (p .- M_d .* d̃) + net[:f̅][i])
        @constraint(model, [μ̲[i], aux_μ̲[i]] in SOS1())
    end
    for i in 1:net[:N]
        @constraint(model, aux_γ̅[i] == net[:p̅][i] - p[i])
        @constraint(model, [γ̅[i], aux_γ̅[i]] in SOS1())
        @constraint(model, aux_γ̲[i] == p[i] - net[:p̲][i])
        @constraint(model, [γ̲[i], aux_γ̲[i]] in SOS1())
    end
    # solve primal-dual model
    optimize!(model)
    @info("done solving the synthetic demand PP problem: $(termination_status(model))")
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        sol = Dict(:status => termination_status(model),
            :obj => JuMP.objective_value(model),
            :p => JuMP.value.(p),
            :opf => net[:c1]' * JuMP.value.(p),
            :CPUtime => solve_time(model),
            :d̃ => JuMP.value.(M_d .* d̃),
            :r_d => JuMP.value.(r_d),
            :r_c => JuMP.value.(r_c),
            :C̃ => C̃_d
        )
    else
        sol = Dict(:status => termination_status(model))
    end
end

# CRO post processing algorithm to generate cyber resilient synthetic load datasets
function CRO_Demand(net,net_c,d̃_0,C̃)
    K = size(net_c[:A],1)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 1)

    # synthetic demand data
    @variable(model, d̃[1:net[:N]])                                      # may not be necessery because DP asnwer can be negative 
    # auxiliary relaxation variable on line capacity
    @variable(model, v[1:net[:E]].>=0)
    # OPF compact parameter matrix
    @variable(model, h[1:K])
    @constraint(model, h[1] == ones(net[:N])'*(M_d.*d̃))
    @constraint(model, h[2:1+net[:E]] .== -net[:F]*(M_d.*d̃) .- net[:f̅] .-v)
    @constraint(model, h[2+net[:E]:K] .== net[:F]*(M_d.*d̃) .- net[:f̅] .-v)
    # power generation
    @variable(model, p1[1:net[:N]])
    @variable(model, p2[1:net[:N]])
    # auxiliary variable for objectives
    @variable(model, r_d)
    @variable(model, r_c_attack)
    @variable(model, r_c_OPF)
    @constraint(model, [r_d ; M_d.*(d̃ .- d̃_0)] in MOI.NormOneCone(1 + length(d̃)))
    @constraint(model, aux_con1, [r_c_attack ; C̃ .- net[:c1]'*p1 .- att[:Ψ]*(ones(net[:E])'*v)] in MOI.NormOneCone(2))
    @constraint(model, aux_con2, [r_c_OPF ; C̃ .- net[:c1]'*p2] in MOI.NormOneCone(2))
    # Cyber Resilient Synthetic Dataset Generation objective function 
    @objective(model, Min, syn_set[:β]*r_c_attack + (1-syn_set[:β])*r_c_OPF + syn_set[:γ]*r_d)

    #-----------------------KKT Conditions of Bilevel Attack Model----------------------
    # primal variables 
    @variable(model, μ̅1[1:K,1:net[:N]].>=0)
    @variable(model, μ̲1[1:K,1:net[:N]].>=0)
    @variable(model, λ1[1:K])
    # primal constraint
    @constraint(model, con_c[k = 1:K], μ̅1[k,:]'*att[:δ̅] .- μ̲1[k,:]'*att[:δ̲]  .<= - net_c[:A][k,:]'*p1 .- h[k])
    @constraint(model, con_d[k = 1:K],  net_c[:B][k,:] .- μ̅1[k,:] .+ μ̲1[k,:] .- ones(net[:N]) .* λ1[k] .== 0)
    @constraint(model, p1 .<= net[:p̅])
    @constraint(model, p1 .>= net[:p̲])

    #dual variables
    @variable(model, θ1[1:K].>=0)
    @variable(model, ζ1[1:K,1:net[:N]])
    @variable(model, ω̅1[1:net[:N]].>=0)
    @variable(model, ω̲1[1:net[:N]].>=0)
    @variable(model, π̅1[1:K,1:net[:N]].>=0)
    @variable(model, π̲1[1:K,1:net[:N]].>=0)
    @variable(model, ξ[1:net[:E]].>=0)

    # stationarity condition
    @constraint(model, con_h, net[:c1] .+ net_c[:A]'θ1 .+ ω̅1 .- ω̲1 .==0)
    @constraint(model, con_i[k=1:K], θ1[k].*att[:δ̅] .-ζ1[k,:] .-π̅1[k,:] .==0)
    @constraint(model, con_j[k=1:K], -θ1[k].*att[:δ̲] .+ζ1[k,:] .-π̲1[k,:] .==0)
    @constraint(model, con_k[k=1:K], ones(net[:N])'*ζ1[k,:]==0)
    @constraint(model, con_l[k=1:K], ones(net[:E]).*att[:Ψ] .- θ1[2:net[:E]+1] .- θ1[net[:E]+2:K] .- ξ .==0)

    # complementarity condition
    @variable(model, aux_θ1[1:K].>=0)
    @variable(model, aux_ω̅1[1:net[:N]].>=0)
    @variable(model, aux_ω̲1[1:net[:N]].>=0)

    @constraint(model, con_aux1[k=1:K] , aux_θ1[k] == -μ̅1[k,:]'*att[:δ̅] + μ̲1[k,:]'*att[:δ̲] - net_c[:A][k,:]'*p1 - h[k])
    @constraint(model, con_m[k=1:K], [θ1[k] , aux_θ1[k]] in SOS1())
    @constraint(model, con_aux2 , aux_ω̅1 .== -p1 .+ net[:p̅])
    @constraint(model, con_n[n=1:net[:N]], [ω̅1[n], aux_ω̅1[n]] in SOS1())
    @constraint(model, con_aux3 , aux_ω̲1 .==  -net[:p̲] .+ p1)
    @constraint(model, con_o[n=1:net[:N]], [ω̲1[n], aux_ω̲1[n]] in SOS1())

    @constraint(model, con_p[k=1:K,n=1:net[:N]], [π̅1[k,n], μ̅1[k,n]] in SOS1())
    @constraint(model, con_q[k=1:K,n=1:net[:N]], [π̲1[k,n], μ̲1[k,n]] in SOS1())
    @constraint(model, con_r[l=1:net[:E]], [ξ[l],v[l]] in SOS1())

    #-----------------------KKT Conditions of Normal OPF----------------------

    # dual variables
    @variable(model, λ2)
    @variable(model, μ̅2[1:net[:E]]>=0)
    @variable(model, μ̲2[1:net[:E]]>=0)
    @variable(model, γ̅2[1:net[:N]]>=0)
    @variable(model, γ̲2[1:net[:N]]>=0)

    # primal constriants 
    @constraint(model, ones(net[:N])'*(p2 .- M_d.*d̃) .== 0)
    @constraint(model, - net[:F]*(p2 .- M_d.*d̃) .>= - net[:f̅])
    @constraint(model, net[:f̅] .>= - net[:F]*(p2 .- M_d.*d̃))
    @constraint(model, -p2 .>= -net[:p̅])
    @constraint(model, -net[:p̲] .>= -p2)

    # stationary constraints
    @constraint(model, net[:c1] + λ2 * ones(net[:N]) + net[:F]' * μ̅2 - net[:F]' * μ̲2 + γ̅2 - γ̲2 .== 0)
    # # strong duality 
    # complementarity conditions formulated as SOS1 constraints
    @variable(model,aux_μ̅2[1:net[:E]].>=0)
    @variable(model,aux_μ̲2[1:net[:E]].>=0)
    @variable(model,aux_γ̅2[1:net[:E]].>=0)
    @variable(model,aux_γ̲2[1:net[:E]].>=0)
    for i in 1:net[:E]
        @constraint(model, aux_μ̅2[i] .== net[:f̅][i] .- net[:F][i,:]'*(p2 .- M_d.*d̃))
        @constraint(model, [μ̅2[i] , aux_μ̅2[i]] in SOS1())

        @constraint(model, aux_μ̲2[i] .== net[:F][i,:]'*(p2 .- M_d.*d̃) .+ net[:f̅][i])
        @constraint(model, [μ̲2[i] , aux_μ̲2[i]] in SOS1())
    end
    for i in 1:net[:N]
        @constraint(model, aux_γ̅2[i] .== net[:p̅][i] .- p2[i])
        @constraint(model, [γ̅2[i], aux_γ̅2[i]] in SOS1())

        @constraint(model, aux_γ̲2[i] .== p2[i] .- net[:p̲][i])
        @constraint(model, [γ̲2[i] , aux_γ̲2[i]] in SOS1())
    end

    # solve primal-dual model
    optimize!(model)
    @info("done solving the Feasible CRData_Demand_Tradeoff problem: $(termination_status(model))")
    status = "$(termination_status(model))"
    if status == "OPTIMAL"
        sol = Dict(:status => termination_status(model),
        :obj => JuMP.objective_value(model),
        :p1 => JuMP.value.(p1),
        :p2 => JuMP.value.(p2),
        :CPUtime => solve_time(model),
        :d̃ => JuMP.value.(M_d.*d̃),
        :r_d => JuMP.value.(r_d),
        :r_c_attack => JuMP.value.(r_c_attack),
        :r_c_OPF => JuMP.value.(r_c_OPF),
        :v => JuMP.value.(v),
        :C̃ => C̃
        )
    else
        sol = Dict(:status => termination_status(model))
    end
end

# CRO post processing algorithm with a specific attack constraint set W_set
function CRO_Exp_Demand(net,net_c,d̃_0,C̃,W_set)
    K = size(net_c[:A],1)
    model = Model(optimizer_with_attributes(Gurobi.Optimizer))
    set_optimizer_attribute(model, "OutputFlag", 1)

    # synthetic demand data
    @variable(model, d̃[1:net[:N]].>=0)                                      # may not be necessery because DP asnwer can be negative 
    # auxiliary relaxation variable on line capacity
    @variable(model, v[1:net[:E]].>=0)
    @variable(model, v2[1:net[:E]].>=0)
    # OPF compact parameter matrix
    @variable(model, h[1:K])
    @constraint(model, h[1] == ones(net[:N])'*(M_d.*d̃))
    @constraint(model, h[2:1+net[:E]] .== -net[:F]*(M_d.*d̃) .- net[:f̅] .-v)
    @constraint(model, h[2+net[:E]:K] .== net[:F]*(M_d.*d̃) .- net[:f̅] .-v)
    # power generation
    @variable(model, p1[1:net[:N]])
    @variable(model, p2[1:net[:N]])
    # auxiliary variable for objectives
    @variable(model, r_d)
    @variable(model, r_c_attack)
    @variable(model, r_c_OPF)
    @constraint(model, [r_d ; M_d.*(d̃ .- d̃_0)] in MOI.NormOneCone(1 + length(d̃)))
    if syn_set[:obj_type]==1
        @constraint(model, aux_con1, [r_c_attack ; C̃ .- net[:c1]'*p1 .- att[:Ψ]*(ones(net[:E])'*v)] in MOI.NormOneCone(2))
    else
        @constraint(model, aux_con1, [r_c_attack ; net[:c1]'*p1 .+ att[:Ψ]*(ones(net[:E])'*v) .- net[:c1]'*p2- att[:Ψ]*(ones(net[:E])'*v2)] in MOI.NormOneCone(2))
    end
    @constraint(model, aux_con2, [r_c_OPF ; C̃ .- net[:c1]'*p2- att[:Ψ]*(ones(net[:E])'*v2)] in MOI.NormOneCone(2))
    # Cyber Resilient Synthetic Dataset Generation objective function 
    @objective(model, Min, syn_set[:β]*r_c_attack + r_c_OPF + syn_set[:γ]*r_d)

    #-----------------------KKT Conditions of Bilevel Attack Model----------------------
    # primal variables 
    @variable(model, μ̅1[1:K,1:net[:N]].>=0)
    @variable(model, μ̲1[1:K,1:net[:N]].>=0)
    @variable(model, λ1[1:K])
    # dual variables
    @variable(model, θ1[1:K].>=0)
    @variable(model, ζ1[1:K,1:net[:N]])
    @variable(model, ω̅1[1:net[:N]].>=0)
    @variable(model, ω̲1[1:net[:N]].>=0)
    @variable(model, π̅1[1:K,1:net[:N]].>=0)
    @variable(model, π̲1[1:K,1:net[:N]].>=0)
    @variable(model, ξ[1:net[:E]].>=0)
    # auxiliary variables for complementarity condition
    @variable(model, aux_θ1[1:K].>=0)
    @variable(model, aux_ω̅1[1:net[:N]].>=0)
    @variable(model, aux_ω̲1[1:net[:N]].>=0)
    for k in 1:K
        if k in W_set
            # primal constraint for attacked constraints
            @constraint(model, μ̅1[k,:]'*att[:δ̅] .- μ̲1[k,:]'*att[:δ̲]  .<= - net_c[:A][k,:]'*p1 .- h[k])
            @constraint(model, net_c[:B][k,:] .- μ̅1[k,:] .+ μ̲1[k,:] .- ones(net[:N]) .* λ1[k] .== 0)
            # stationarity condition for attacked constraints
            @constraint(model, net[:c1] .+ net_c[:A]'θ1 .+ ω̅1 .- ω̲1 .==0)
            @constraint(model, θ1[k].*att[:δ̅] .-ζ1[k,:] .-π̅1[k,:] .==0)
            @constraint(model, -θ1[k].*att[:δ̲] .+ζ1[k,:] .-π̲1[k,:] .==0)
            @constraint(model, ones(net[:N])'*ζ1[k,:]==0)
            # complementarity condition for attacked constraints
            @constraint(model, aux_θ1[k] == -μ̅1[k,:]'*att[:δ̅] + μ̲1[k,:]'*att[:δ̲] - net_c[:A][k,:]'*p1 - h[k])
            @constraint(model,  [θ1[k] , aux_θ1[k]] in SOS1())

            for n in 1:net[:N]
                @constraint(model,  [π̅1[k,n], μ̅1[k,n]] in SOS1())
                @constraint(model,  [π̲1[k,n], μ̲1[k,n]] in SOS1())
            end
        else
            @constraint(model, 0 .<= - net_c[:A][k,:]'*p1 .- h[k])
        end
    end
    # common primal constraint
    @constraint(model, p1 .<= net[:p̅])
    @constraint(model, p1 .>= net[:p̲])
    # common stationarity condition
    @constraint(model, ones(net[:E]).*att[:Ψ] .- θ1[2:net[:E]+1] .- θ1[net[:E]+2:K] .- ξ .==0)
    # common complimentary constraints
    @constraint(model, con_aux2 , aux_ω̅1 .== -p1 .+ net[:p̅])
    @constraint(model, con_n[n=1:net[:N]], [ω̅1[n], aux_ω̅1[n]] in SOS1())
    @constraint(model, con_aux3 , aux_ω̲1 .==  -net[:p̲] .+ p1)
    @constraint(model, con_o[n=1:net[:N]], [ω̲1[n], aux_ω̲1[n]] in SOS1())
    @constraint(model, con_r[l=1:net[:E]], [ξ[l],v[l]] in SOS1())
            
    #-----------------------KKT Conditions of Normal OPF----------------------

    # dual variables
    @variable(model, λ2)
    @variable(model, μ̅2[1:net[:E]]>=0)
    @variable(model, μ̲2[1:net[:E]]>=0)
    @variable(model, γ̅2[1:net[:N]]>=0)
    @variable(model, γ̲2[1:net[:N]]>=0)
    @variable(model, ν[1:net[:E]]>=0) # dual for v>=0


    # primal constriants 
    @constraint(model, ones(net[:N])'*(p2 .- M_d.*d̃) .== 0)
    @constraint(model, - net[:F]*(p2 .- M_d.*d̃) .>= - net[:f̅].-v2)
    @constraint(model, net[:f̅].+v2 .>= - net[:F]*(p2 .- M_d.*d̃))
    @constraint(model, -p2 .>= -net[:p̅])
    @constraint(model, -net[:p̲] .>= -p2)

    # stationary constraints
    @constraint(model, net[:c1] + λ2 * ones(net[:N]) + net[:F]' * μ̅2 - net[:F]' * μ̲2 + γ̅2 - γ̲2 .== 0)
    @constraint(model, att[:Ψ]*(ones(net[:E])) .- μ̅2 .- μ̲2 .- ν== 0)
    # # strong duality 
    # complementarity conditions formulated as SOS1 constraints
    @variable(model,aux_μ̅2[1:net[:E]].>=0)
    @variable(model,aux_μ̲2[1:net[:E]].>=0)
    @variable(model,aux_γ̅2[1:net[:E]].>=0)
    @variable(model,aux_γ̲2[1:net[:E]].>=0)
    @variable(model,aux_ν[1:net[:E]].>=0)
    for i in 1:net[:E]
        @constraint(model, aux_μ̅2[i] .== net[:f̅][i].+v2[i] .- net[:F][i,:]'*(p2 .- M_d.*d̃))
        @constraint(model, [μ̅2[i] , aux_μ̅2[i]] in SOS1())

        @constraint(model, aux_μ̲2[i] .== net[:F][i,:]'*(p2 .- M_d.*d̃) .+ net[:f̅][i].+v2[i])
        @constraint(model, [μ̲2[i] , aux_μ̲2[i]] in SOS1())

        @constraint(model, aux_ν[i] .== v2[i])
        @constraint(model, [ν[i] , aux_ν[i]] in SOS1())
    end
    for i in 1:net[:N]
        @constraint(model, aux_γ̅2[i] .== net[:p̅][i] .- p2[i])
        @constraint(model, [γ̅2[i], aux_γ̅2[i]] in SOS1())

        @constraint(model, aux_γ̲2[i] .== p2[i] .- net[:p̲][i])
        @constraint(model, [γ̲2[i] , aux_γ̲2[i]] in SOS1())
    end

    # solve primal-dual model
    optimize!(model)
    @info("done solving the Feasible CRData_Demand_Tradeoff problem: $(termination_status(model))")
    status = "$(termination_status(model))"
    if status == "OPTIMAL" || status == "LOCALLY_SOLVED"
        sol = Dict(:status => termination_status(model),
        :obj => JuMP.objective_value(model),
        :p1 => JuMP.value.(p1),
        :p2 => JuMP.value.(p2),
        :CPUtime => solve_time(model),
        :d̃ => JuMP.value.(M_d.*d̃),
        :r_d => JuMP.value.(r_d),
        :r_c_attack => JuMP.value.(r_c_attack),
        :r_c_OPF => JuMP.value.(r_c_OPF),
        :v => JuMP.value.(v),
        :C̃ => C̃
        )
    else
        sol = Dict(:status => termination_status(model))
    end
end

# Use exponential mechanism to select a number of con_num worst-cast constraint under attack
function Exponential_Mechanism(net,net_c,con_num=10)
    # Random.seed!(33)
    K = size(net_c[:A],1)
    unattacked_set=collect(1:K)
    attacked_list=[]
    ΔC_max_list=[]
    net_test=deepcopy(net_c)
    net_test[:c]=vcat(ones(net[:N])'*net[:d] , -net[:F]*net[:d]-net[:f̅] , net[:F]*net[:d]-net[:f̅])
    sol_opf_k=DC_OPF_compact_notation(net_c)
    while length(unattacked_set) != K-con_num
        # find the first worst case attack scenarios
        if length(unattacked_set)==K
            ΔC_w=zeros(K)
            for i in 1:K
                sol_attack_k=AttackModel_RO_k(net,net_test,att,unattacked_set[i])
                ΔC_w[i]=sol_attack_k[:obj]+rand(Laplace(0,(syn_set[:c_max]*syn_set[:α])/(syn_set[:ϵ]/3/att[:τ])))
                # ΔC_w[i]=abs(sol_attack_k[:obj]-sol_opf_k[:obj])+rand(Laplace(0,(syn_set[:c_max]*syn_set[:α])/(syn_set[:ϵ]/3/att[:τ])))
            end
            ΔC_max , w_max_idx=findmax(ΔC_w)
            push!(ΔC_max_list,ΔC_max)
            # add the worst case attack in this round to the attacked list
            push!(attacked_list,unattacked_set[w_max_idx])
            # delete the constraint that's been added to attacked list
            deleteat!(unattacked_set,w_max_idx)
        else
            unattacked_constraint_num=length(unattacked_set)
            ΔC_w=zeros(unattacked_constraint_num)
            for i in 1:unattacked_constraint_num
                sol_attack_k=AttackModel_RO_k(net,net_test,att,vcat(attacked_list,unattacked_set[i]))
                # ΔC_w[i]=abs.(sol_attack_k[:obj]-sol_opf_k[:obj])+rand(Laplace(0,(syn_set[:c_max]*syn_set[:α])/syn_set[:ϵ]))
                ΔC_w[i]=sol_attack_k[:obj]+rand(Laplace(0,(syn_set[:c_max]*syn_set[:α])/syn_set[:ϵ]/3/att[:τ]))
            end
            ΔC_max , w_max_idx=findmax(ΔC_w)
            push!(ΔC_max_list,ΔC_max)
            push!(attacked_list,unattacked_set[w_max_idx])
            deleteat!(unattacked_set,w_max_idx)
        end 
    end
    return attacked_list
end

# case study of CRO Algorithms in distributions
function Results_CRO(net,net_c,sample_points=10)
    # calculate the standard OPF cost
    sol_opf=DC_OPF_compact_feas(net_c)
    # println("Standard OPF Objectives: ",sol_opf[:obj]) 
    dict_pp=Dict(:C_OPF => zeros(sample_points), :C_att_RO => zeros(sample_points), :C_att_RO_τ1 => zeros(sample_points),:C_att_RO_τ5 => zeros(sample_points),:C_att_BO => zeros(sample_points), :d̃ => zeros(sample_points,net[:N]), :d̃0 => zeros(sample_points,net[:N]),:v => zeros(sample_points,net[:E]), :C̃ =>zeros(sample_points))
    dict_tradeoff=Dict(:C_OPF => zeros(sample_points), :C_att_RO => zeros(sample_points), :C_att_RO_τ1 => zeros(sample_points),:C_att_RO_τ5 => zeros(sample_points),:C_att_BO => zeros(sample_points), :d̃ => zeros(sample_points,net[:N]),:d̃0 => zeros(sample_points,net[:N]),:v => zeros(sample_points,net[:E]), :C̃ =>zeros(sample_points))
    Random.seed!(30)
    net_test=deepcopy(net)
    net_test[:A]=net_c[:A]
    net_test[:B]=net_c[:B]
    # @showprogress for i in 1:sample_points
    for i in 1:sample_points
        d̃_i= M_d.*(net[:d] .+ rand(Laplace(0,syn_set[:α]/(syn_set[:ϵ]/2)),net[:N]))
        C̃_i=sol_opf[:obj] .+ rand(Laplace(0,(syn_set[:c_max]*syn_set[:α])/(syn_set[:ϵ]/2)))
        println("synthetic cost:", C̃_i)
        println("synthetic demand:", d̃_i)
        # generate dataset using standard post processing
        sol_pp_i=OPF_PP_Demand(net,d̃_i,C̃_i)
        net_test[:c]=vcat(ones(net[:N])'*sol_pp_i[:d̃] , -net[:F]*sol_pp_i[:d̃]-net[:f̅] , net[:F]*sol_pp_i[:d̃]-net[:f̅])
        # Record the OPF objectives on the generated synthetic datasets
        sol_opf_pp_i=DC_OPF_compact_feas(net_test)
        dict_pp[:C_OPF][i]=sol_opf_pp_i[:obj]
        dict_pp[:v][i,:]=sol_opf_pp_i[:v] 
        # Record the synthetic demand data using standard postporocessing
        dict_pp[:d̃][i,:]=sol_pp_i[:d̃]
        # Record the original demand and OPF cost after obfuscation
        dict_pp[:d̃0][i,:]=d̃_i
        dict_pp[:C̃][i]=C̃_i
        # record the Attack objectives on the generated synthetic datasets
        sol_att_pp_i=AttackModel_RO(net,net_test,att)
        dict_pp[:C_att_RO][i]=sol_att_pp_i[:obj]
        net_test[:d]=sol_pp_i[:d̃]
        dict_pp[:C_att_BO][i]=AttackModel_BO_feas(net_test,att)[:obj]

        # att damage of τ-constraint attacks 
        # attack_list_τ=Exponential_Mechanism(net,net_c,5)
        # dict_pp[:C_att_RO_τ1][i]=AttackModel_feas_k(net, net_test,att,attack_list_τ[1])[:obj]
        # dict_pp[:C_att_RO_τ5][i]=AttackModel_feas_k(net, net_test,att,attack_list_τ)[:obj]

        # Record the synthetic demand data using tradeoff postporocessing model
        sol_tradeoff_i=CRO_Demand(net,net_c,d̃_i,C̃_i)
        dict_tradeoff[:d̃][i,:]=sol_tradeoff_i[:d̃]
        # Record the OPF objectives on the generated synthetic datasets
        net_test[:c]=vcat(ones(net[:N])'*sol_tradeoff_i[:d̃] , -net[:F]*sol_tradeoff_i[:d̃]-net[:f̅] , net[:F]*sol_tradeoff_i[:d̃]-net[:f̅])
        sol_opf_tradeoff_i=DC_OPF_compact_feas(net_test)
        dict_tradeoff[:C_OPF][i]=sol_opf_tradeoff_i[:obj]
        dict_tradeoff[:v][i,:]=sol_opf_tradeoff_i[:v]  
        dict_tradeoff[:d̃0][i,:]=d̃_i
        dict_tradeoff[:C̃][i]=C̃_i
        # record the attack objectives on the generated synthetic datasets
        dict_tradeoff[:C_att_RO][i]=AttackModel_RO(net,net_test,att)[:obj]
        net_test[:d]=sol_tradeoff_i[:d̃]
        dict_tradeoff[:C_att_BO][i]=AttackModel_BO_feas(net_test,att)[:obj]

        # dict_tradeoff[:C_att_RO_τ1][i]=AttackModel_feas_k(net, net_test,att,attack_list_τ[1])[:obj]
        # dict_tradeoff[:C_att_RO_τ5][i]=AttackModel_feas_k(net, net_test,att,attack_list_τ)[:obj]
    end
    return dict_pp, dict_tradeoff
end

# serialize the results of CRO and save it results_CRO folder
function CaseStudy_CRO(sample_num=10)
    syn_set[:α]= 50
    dict_pp_50 , dict_cro_50 = Results_CRO(sample_num)
    @show [mean(dict_cro_50[:C_OPF]) std(dict_cro_50[:C_OPF]) mean(dict_cro_50[:C_att_BO]) std(dict_cro_50[:C_att_BO]) mean(dict_pp_50[:C_OPF]) std(dict_pp_50[:C_OPF]) mean(dict_pp_50[:C_att_BO]) std(dict_pp_50[:C_att_BO])]
    # save the results to file
    serialize("./results_CRO_feas/dict_pp_50.jls", dict_pp_50)
    serialize("./results_CRO_feas/dict_cro_50.jls", dict_cro_50)

    syn_set[:α]= 100
    dict_pp_100 , dict_cro_100 = Results_CRO(sample_num)
    @show [mean(dict_cro_100[:C_OPF]) std(dict_cro_100[:C_OPF]) mean(dict_cro_100[:C_att_BO]) std(dict_cro_100[:C_att_BO]) mean(dict_pp_100[:C_OPF]) std(dict_pp_100[:C_OPF]) mean(dict_pp_100[:C_att_BO]) std(dict_pp_100[:C_att_BO])]
    # save the results to file
    serialize("./results_CRO_feas/dict_pp_100.jls", dict_pp_100)
    serialize("./results_CRO_feas/dict_cro_100.jls", dict_cro_100)

    syn_set[:α]= 20
    dict_pp_20 , dict_cro_20 = Results_CRO(sample_num)
    @show [mean(dict_cro_20[:C_OPF]) std(dict_cro_20[:C_OPF]) mean(dict_cro_20[:C_att_BO]) std(dict_cro_20[:C_att_BO]) mean(dict_pp_20[:C_OPF]) std(dict_pp_20[:C_OPF]) mean(dict_pp_20[:C_att_BO]) std(dict_pp_20[:C_att_BO])]
    # save the results to file
    serialize("./results_CRO_feas/dict_pp_20.jls", dict_pp_20)
    serialize("./results_CRO_feas/dict_cro_20.jls", dict_cro_20)
end

# figures of CRO algorithms
function plot_CRO_test(net,net_c,sample_num=50)
    # draw a distribution of opf objectives before and after attack 
    # Create a figure
    syn_set[:α]=20
    dict_pp_20,dict_cro_20=Results_CRO(net,net_c,sample_num)
    syn_set[:α]=50
    dict_pp_50,dict_cro_50=Results_CRO(net,net_c,sample_num)
    syn_set[:α]=100
    dict_pp_100,dict_cro_100=Results_CRO(net,net_c,sample_num)

    fig = Figure(size=(1080, 750))
    transparency=0.4
    bin_num=10

    ax1 = Axis(fig[1, 1], title="Adjacency α=20MW", xticks=xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false,xticksvisible = false,xticklabelsvisible = false,ylabel="Frequency",titlesize=29, ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    # Plot subfugure (1,1)
    hist!(ax1, dict_pp_20[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{opf}(\tilde{\textbf{d}}_\text{pp})")
    hist!(ax1, dict_pp_20[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{pp})")
    vlines!(ax1, mean(dict_pp_20[:C_OPF]), linestyle=:dot, linewidth=2, color="#003049")
    vlines!(ax1, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax1, mean(dict_pp_20[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")


    ax2 = Axis(fig[2, 1],xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false,xlabel="Dispatch Cost(k\$)",ylabel="Frequency",titlesize=29,ylabelsize=24,xlabelsize=24, xticklabelsize=20, aspect=1)
    hist!(ax2, dict_cro_20[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{opf}(\tilde{\textbf{d}}_\text{cro})")
    hist!(ax2, dict_cro_20[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{cro})")
    vlines!(ax2, mean(dict_cro_20[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax2, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax2, mean(dict_cro_20[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")

    ax3 = Axis(fig[1, 2], title="Adjacency α=50MW",xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)), xticksvisible = false,xticklabelsvisible = false,yticksvisible = false,yticklabelsvisible = false,titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    hist!(ax3, dict_pp_50[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{opf}(\tilde{\textbf{d}}_\text{pp})")
    hist!(ax3, dict_pp_50[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{pp})")
    vlines!(ax3, mean(dict_pp_50[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax3, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax3, mean(dict_pp_50[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")


    ax4 = Axis(fig[2, 2],xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false, xlabel="Dispatch Cost(k\$)", titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    hist!(ax4, dict_cro_50[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{cro})")
    hist!(ax4, dict_cro_50[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{cro})")
    vlines!(ax4, mean(dict_cro_50[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax4, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax4, mean(dict_cro_50[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")

    ax5 = Axis(fig[1, 3], title="Adjacency α=100MW",xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)), xticksvisible = false,xticklabelsvisible = false,titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20,yticksvisible = false,yticklabelsvisible = false, aspect=1)
    hist!(ax5, dict_pp_100[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency))
    hist!(ax5, dict_pp_100[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency))
    vlines!(ax5, mean(dict_pp_100[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax5, [88200], linestyle=:dot, linewidth=3, color="#fca311")
    vlines!(ax5, mean(dict_pp_100[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")
    

    ax6 = Axis(fig[2, 3],xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false, xlabel="Dispatch Cost(k\$)", titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    hist!(ax6, dict_cro_100[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency))
    hist!(ax6, dict_cro_100[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency))
    vlines!(ax6, mean(dict_cro_100[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax6, [88200], linestyle=:dot, linewidth=3, color="#fca311")
    vlines!(ax6, mean(dict_cro_100[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")
    
    axislegend(ax3,labelsize=25,position=:lt)
    axislegend(ax4,labelsize=25,position=:lt)
    # Show the figure
    display(fig)
    # save("CRO_alpha.png",fig)
end

# deserialize the results of CRO for plotting
function Plot_CRO()
    bin_num=35
    transparency=0.4

    # read the file 

    dict_pp_50 = deserialize("./results_CRO/dict_pp_50.jls")
    dict_cro_50 = deserialize("./results_CRO/dict_cro_50.jls")
    dict_pp_100 = deserialize("./results_CRO/dict_pp_100.jls")
    dict_cro_100 = deserialize("./results_CRO/dict_cro_100.jls")
    dict_pp_20 = deserialize("./results_CRO/dict_pp_20.jls")
    dict_cro_20 = deserialize("./results_CRO/dict_cro_20.jls")

    # Create a figure
    fig = Figure(size=(1080, 750))


    ax1 = Axis(fig[1, 1], title="Adjacency α=20MW", xticks=xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false,xticksvisible = false,xticklabelsvisible = false,ylabel="Frequency",titlesize=29, ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    # Plot subfugure (1,1)
    hist!(ax1, dict_pp_20[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{opf}(\tilde{\textbf{d}}_\text{pp})")
    hist!(ax1, dict_pp_20[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{pp})")
    vlines!(ax1, mean(dict_pp_20[:C_OPF]), linestyle=:dot, linewidth=2, color="#003049")
    vlines!(ax1, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax1, mean(dict_pp_20[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")


    ax2 = Axis(fig[2, 1],xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false,xlabel="Dispatch Cost(k\$)",ylabel="Frequency",titlesize=29,ylabelsize=24,xlabelsize=24, xticklabelsize=20, aspect=1)
    hist!(ax2, dict_cro_20[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{opf}(\tilde{\textbf{d}}_\text{cro})")
    hist!(ax2, dict_cro_20[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{cro})")
    vlines!(ax2, mean(dict_cro_20[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax2, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax2, mean(dict_cro_20[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")

    ax3 = Axis(fig[1, 2], title="Adjacency α=50MW",xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)), xticksvisible = false,xticklabelsvisible = false,yticksvisible = false,yticklabelsvisible = false,titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    hist!(ax3, dict_pp_50[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{opf}(\tilde{\textbf{d}}_\text{pp})")
    hist!(ax3, dict_pp_50[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{pp})")
    vlines!(ax3, mean(dict_pp_50[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax3, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax3, mean(dict_pp_50[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")


    ax4 = Axis(fig[2, 2],xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false, xlabel="Dispatch Cost(k\$)", titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    hist!(ax4, dict_cro_50[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{cro})")
    hist!(ax4, dict_cro_50[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency), label=L"C_\text{att}^{\text{BO}}(\tilde{\textbf{d}}_\text{cro})")
    vlines!(ax4, mean(dict_cro_50[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax4, [88200], linestyle=:dot, linewidth=3, color="#fca311",label=L"C_\text{opf}(\textbf{d})")
    vlines!(ax4, mean(dict_cro_50[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")

    ax5 = Axis(fig[1, 3], title="Adjacency α=100MW",xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)), xticksvisible = false,xticklabelsvisible = false,titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20,yticksvisible = false,yticklabelsvisible = false, aspect=1)
    hist!(ax5, dict_pp_100[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency))
    hist!(ax5, dict_pp_100[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency))
    vlines!(ax5, mean(dict_pp_100[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax5, [88200], linestyle=:dot, linewidth=3, color="#fca311")
    vlines!(ax5, mean(dict_pp_100[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")
    

    ax6 = Axis(fig[2, 3],xticks=([75000,80000,85000,90000,95000,100000,105000],["75","80","85","90","95","100","105"]),limits=((70000, 108000),(0,0.00030)),yticksvisible = false,yticklabelsvisible = false, xlabel="Dispatch Cost(k\$)", titlesize=29,ylabelsize=24,xlabelsize=24,xticklabelsize=20, aspect=1)
    hist!(ax6, dict_cro_100[:C_OPF], bins=bin_num, normalization=:pdf, color=("#003049",transparency))
    hist!(ax6, dict_cro_100[:C_att_BO], bins=bin_num, normalization=:pdf, color=("#780000",transparency))
    vlines!(ax6, mean(dict_cro_100[:C_OPF]), linestyle=:dot, linewidth=3, color="#003049")
    vlines!(ax6, [88200], linestyle=:dot, linewidth=3, color="#fca311")
    vlines!(ax6, mean(dict_cro_100[:C_att_BO]), linestyle=:dot, linewidth=3, color="#780000")
    
    axislegend(ax3,labelsize=25,position=:lt)
    axislegend(ax4,labelsize=25,position=:lt)
    # Show the figure
    display(fig)
    save("CRO_alpha.png",fig)
end

# case study of CRO-Exp algorithm, return the attack damage in absolute values
function CaseStudy_CRO_EM(net,net_c,sample_points=100)
    # This functions include all results from τ=1 to τ=att[:τ]
    dict_cr=Dict(:C_OPF => zeros(att[:τ]+1,sample_points), :C_att_RO => zeros(att[:τ]+1,sample_points),:C_att_BO => zeros(att[:τ]+1,sample_points), :d̃ => zeros(sample_points,net[:N]), :d̃0 => zeros(sample_points,net[:N]),:v => zeros(sample_points,net[:E]), :C̃ =>zeros(sample_points))
    Random.seed!(28)
    net_c_test=deepcopy(net_c)
    net_c_test[:c]=vcat(ones(net[:N])'*net[:d] , -net[:F]*net[:d]-net[:f̅] , net[:F]*net[:d]-net[:f̅])
    sol_opf=DC_OPF_compact_feas(net_c)
    net_temp=deepcopy(net)
    
    # Constraint_Set=AttackModel_worst_k(net,net_c,att,att[:τ])
    for i in 1:sample_points
        d̃_i= M_d.*(net[:d] .+ rand(Laplace(0,syn_set[:α]/(syn_set[:ϵ]/3)),net[:N]))
        C̃_i=sol_opf[:obj] .+ rand(Laplace(0,(syn_set[:c_max]*syn_set[:α])/(syn_set[:ϵ]/3)))
        # Record the original demand and OPF cost after obfuscation
        dict_cr[:d̃0][i,:]=d̃_i
        dict_cr[:C̃][i]=C̃_i
        # ensure that test case use the same attack set to show the patterns of τ
        attacked_list_total=Exponential_Mechanism(net,net_c,att[:τ])
        for τ_idx in 0:att[:τ]
            if τ_idx==0
                sol_pp_i=OPF_PP_Demand(net,d̃_i,C̃_i)
                sol_opf_pp_i=OPF(net,sol_pp_i[:d̃])
                dict_cr[:C_OPF][τ_idx+1,i]=sol_opf_pp_i[:obj]
                net_temp[:d]=sol_pp_i[:d̃]
                sol_pp_att_BO=AttackModel_BO(net_temp,att)
                dict_cr[:C_att_BO][τ_idx+1,i]=sol_pp_att_BO[:obj]
            else
                att_list_τ=attacked_list_total[1:τ_idx]
                sol_cr_i=CRO_EM_Demand(net,net_c,d̃_i,C̃_i,att_list_τ)
                # record the Attack objectives on the generated synthetic datasets
                net_c_test[:c]=vcat(ones(net[:N])'*sol_cr_i[:d̃] , -net[:F]*sol_cr_i[:d̃]-net[:f̅] , net[:F]*sol_cr_i[:d̃]-net[:f̅])
                # Record the OPF objectives on the generated synthetic datasets
                sol_opf_cr_i=DC_OPF_compact_feas(net_c_test)
                dict_cr[:C_OPF][τ_idx+1,i]=sol_opf_cr_i[:obj]
                # Record the cost of BO attacks
                net_temp[:d]=sol_cr_i[:d̃]
                sol_cr_att_BO=AttackModel_BO(net_temp,att)
                dict_cr[:C_att_BO][τ_idx+1,i]=sol_cr_att_BO[:obj]
            end
        end
    end
    return dict_cr
end

# case study of CRO-Exp algorithm, return the attack damage in percentages.
function Results_CRO_EM(sample_points=100)
    Random.seed!(8)
    dict_cro_τ=CaseStudy_CRO_EM(net,net_c,sample_points)
    BO_mean_list=[]
    BO_95_list=[]
    BO_05_list=[]
    for τ_idx=1:att[:τ]+1
        println("---------τ=$(τ_idx-1)-----------")
        dict_perc_att_BO=(dict_cro_τ[:C_att_BO][τ_idx,:].-dict_cro_τ[:C_OPF][τ_idx,:])./dict_cro_τ[:C_OPF][τ_idx,:]*100
        dict_perc_att_BO=round.(dict_perc_att_BO,digits=3)
        push!(BO_mean_list,mean(dict_perc_att_BO))
        push!(BO_95_list,quantile(dict_perc_att_BO, 0.9))
        push!(BO_05_list,quantile(dict_perc_att_BO, 0.1))
        @show [mean(dict_perc_att_BO) quantile(dict_perc_att_BO, 0.9) quantile(dict_perc_att_BO, 0.1)]
    end
    println("---------Results-----------")
    @show BO_mean_list
    @show BO_95_list
    @show BO_05_list
end
