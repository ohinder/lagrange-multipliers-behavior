println("loading libaries and functions ...")

using JuMP, MAT, Ipopt, MathProgBase, NLPModels, CUTEst, DataFrames, CSV
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

function add_solver_results!(hist::Array{OnePhase.generic_alg_history,1}, nlp::AbstractNLPModel, inner, t::Int64)
    # TODO
    # check this carefully

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
        if lcon[i] != ucon[i] && !isinf(lcon[i]) && !isinf(ucon[i])
            warn("Please. Write split inequalities so strict complementarity can be correctly computed. JuMP interface currently makes this difficult to do unless inequalities are split.")
        end
    end

    r_U = ucon - g_val
    r_L = g_val - lcon
    # TODO check carefully -- only seems to work with equality constraints with r.h.s of zero????
    con_vio = max(maximum(-r_U),maximum(-r_L),maximum(-s_L),maximum(-s_U),0.0)

    # doesn't allow inequality constraints.
    comp_vec = [mult_x_L .* max.(0.0,x-nlp.meta.lvar); mult_x_U .* max.(0.0,nlp.meta.uvar-x)]
    comp_vec[[nlp.meta.lvar; nlp.meta.uvar] .== Inf] = 0.0
    comp_vec[[nlp.meta.lvar; nlp.meta.uvar] .== -Inf] = 0.0
    #@show comp_vec
    comp = maximum(comp_vec)
    max_comp = maximum(comp_vec)
    min_comp = minimum(comp_vec)

    fval = obj(nlp,x)
    strict_comp_r_U = (r_U + mult_g)[ucon .!= lcon]
    strict_comp_r_L = (r_L + mult_g)[ucon .!= lcon]
    strict_comp = min(minimum(mult_x_L + s_L),minimum(mult_x_U + s_U))

    this_it = OnePhase.generic_alg_history(t,fval,norm_grad_lag,comp,con_vio,con_vio,y_norm,x_norm,strict_comp)

    push!(hist,this_it)
end

function DictToKeyValue(dic::Dict)
    return [(val,key) for (val,key) in dic]
end

function IPOPT_solver_history(model_builder::Function, solver_info::Dict; print_level=0::Int64)
    temp_m = model_builder()
    nlp = NLPModels.MathProgNLPModel(temp_m)
    status = NaN;

    hist = Array{OnePhase.generic_alg_history,1}()
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
    hist = Array{OnePhase.generic_alg_history,1}()
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

function Record_solver_histories(solver_dic::Dict, build_nlp::Function)
    hist_dic = Dict{String,Array{OnePhase.abstract_alg_history,1}}()
    status_dic = Dict{String,Symbol}()
    for (solver_name,solver_info) in solver_dic
        if solver_info["solver"] == :Ipopt
            hist_dic[solver_name], status_dic[solver_name] = IPOPT_solver_history(build_nlp, solver_info);

        elseif solver_info["solver"] == :OnePhase
            m = build_nlp()
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
