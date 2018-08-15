println("loading libaries and functions ...")

using JuMP, MAT, Ipopt, MathProgBase, NLPModels, CUTEst, DataFrames, CSV, PyPlot
using Gurobi # to check if LP is feasible or not
using Suppressor # to stop IPOPT showing a million warnings

if ENV["USER"] == "Oliver" # me
    # so I can work with the latest one-phase IPM
    include("../../../one-phase-2.0/src/OnePhase.jl")
else
    # everyone else: download it from https://github.com/ohinder/OnePhase.jl
    using OnePhase
end
include("lp.jl")
include("plot.jl")

function build_solver(solver_info::Dict)
    if solver_info["solver"] == :OnePhase
        return OnePhase.OnePhaseSolver(solver_info["options"])
    elseif solver_info["solver"] == :Ipopt
        return IpoptSolver(solver_info["options"])
    else
        error("Unknown solver.")
    end
end

struct ipopt_alg_history <: OnePhase.abstract_alg_history
    t::Int64
    #mu::Float64
    fval::Float64
    #dual_scaled::Float64
    norm_grad_lag::Float64
    sy_inf::Float64
    #comp_ratio::Float64
    #sy_mean::Float64
    #primal_residual::Float64
    #farkas::Float64
    #delta::Float64
    #eval_merit_function::Float64
    #eval_phi::Float64
    #eval_grad_phi::Float64
    primal_residual::Float64
    con_vio::Float64
    y_norm::Float64
    x_norm::Float64
    strict_comp::Float64
end

function add_solver_results!(hist::Array{ipopt_alg_history,1}, nlp::AbstractNLPModel, inner, t::Int64)
    # this is where we compute primal feasibility, dual feasibility, etc of Ipopt.
    x = inner.x
    mult_x_L = inner.mult_x_L
    mult_x_U = inner.mult_x_U
    mult_g = inner.mult_g

    x_norm = norm(x,Inf)

    y_norm = max(norm(mult_g,Inf), norm(mult_x_L,Inf), norm(mult_x_U,Inf))
    J = jac(nlp, x)
    norm_grad_lag = norm(grad(nlp,x) +  J' *  mult_g - mult_x_L + mult_x_U,Inf)

    g_val = cons(nlp, x)
    g_val = inner.g

    s_U = nlp.meta.uvar - x
    s_L = x - nlp.meta.lvar

    lcon = nlp.meta.lcon
    ucon = nlp.meta.ucon

    for i = 1:length(lcon)
        if lcon[i] != ucon[i]
            warn("Please turn inequalities into equalities using slack variables. JuMP interface currently makes this difficult handle otherwise and this is essentially what Ipopt does internally.")
        end
    end

    r_U = ucon - g_val
    r_L = g_val - lcon
    con_vio = max(maximum(-r_U),maximum(-r_L),maximum(-s_L),maximum(-s_U),0.0)

    sy_vec = [mult_x_L .* max.(0.0,s_L); mult_x_U .* max.(0.0,s_U)]
    sy_vec[[nlp.meta.lvar; nlp.meta.uvar] .== Inf] = 0.0
    sy_vec[[nlp.meta.lvar; nlp.meta.uvar] .== -Inf] = 0.0
    sy_inf = maximum(sy_vec)

    fval = obj(nlp,x)

    strict_comp = min(minimum(mult_x_L + s_L),minimum(mult_x_U + s_U))

    this_it = ipopt_alg_history(t,fval,norm_grad_lag,sy_inf,con_vio,con_vio,y_norm,x_norm,strict_comp)

    push!(hist,this_it)
end

function DictToKeyValue(dic::Dict)
    return [(val,key) for (val,key) in dic]
end

function IPOPT_solver_history(model_builder::Function, solver_info::Dict; print_level=0::Int64)
    temp_m = model_builder()
    nlp = NLPModels.MathProgNLPModel(temp_m)
    status = NaN;

    hist = Array{ipopt_alg_history,1}()
    opts = solver_info["options"]
    for k = 0:opts[:max_iter]
        m = model_builder()
        new_opts = deepcopy(opts)
        new_opts[:max_iter] = k
        new_opts[:print_level] = print_level
        tmp_solver = IpoptSolver(DictToKeyValue(new_opts))
        setsolver(m, tmp_solver)
        status = @suppress_err solve(m)
        inner = m.internalModel.inner
        add_solver_results!(hist, nlp, inner, k)
        if status != :UserLimit
            break
        end
    end

    if opts[:print_level] > print_level
        # run IPOPT as normal and print output
        m = model_builder()
        tmp_solver = IpoptSolver(DictToKeyValue(opts))
        setsolver(m, tmp_solver)
        status = solve(m)
    end

    return hist,status
end

function IPOPT_solver_history(nlp::NLPModels.AbstractNLPModel, solver_info::Dict; print_level=0::Int64) #max_it::Int64; bound_relax_factor::Float64=0.0,tol::Float64=1e-8)
    hist = Array{ipopt_alg_history,1}()
    opts = solver_info["options"]
    for k = 0:opts[:max_iter]
        if k % mod == 0
            new_opts = deepcopy(opts)
            new_opts[:max_iter] = k
            new_opts[:print_level] = print_level
            tmp_solver = IpoptSolver(DictToKeyValue(new_opts))
            mp = NLPtoMPB(nlp,tmp_solver)
            status = @suppress_err MathProgBase.optimize!(mp)
            inner = mp.inner
            add_solver_results!(hist, nlp, inner, k)
            if status != -1
                break
            end
        end
    end

    return hist
end

function Record_solver_histories(solver_dic::Dict, build_nlp_OnePhase::Function, build_nlp_Ipopt::Function)
    hist_dic = Dict{String,Array{OnePhase.abstract_alg_history,1}}()
    status_dic = Dict{String,Symbol}()
    for (solver_name,solver_info) in solver_dic
        if solver_info["solver"] == :Ipopt
            hist_dic[solver_name], status_dic[solver_name] = IPOPT_solver_history(build_nlp_Ipopt, solver_info);

        elseif solver_info["solver"] == :OnePhase
            m = build_nlp_OnePhase()
            solver = build_solver(solver_info)
            setsolver(m,solver)
            status_dic[solver_name] = solve(m)
            hist_dic[solver_name] = OnePhase.major_its_only(m.internalModel.inner.hist);
        else
            error("Unknown solver.")
        end
    end
    return hist_dic, status_dic
end
