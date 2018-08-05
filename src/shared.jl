using JuMP, MAT, Ipopt, MathProgBase, NLPModels, CUTEst
using Gurobi # to check if LP is feasible or not
#using OnePhase
include("../../../one-phase-2.0/src/OnePhase.jl")
include("lp.jl")
include("plot.jl")

function add_solver_results!(hist::Array{OnePhase.generic_alg_history,1}, nlp::AbstractNLPModel, inner, t::Int64)
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
    # only seems to work with equality constraints with r.h.s of zero.
    con_vio = max(maximum(g_val-nlp.meta.ucon),maximum(nlp.meta.lcon-g_val),0.0)
    #@show nlp.meta.ucon, nlp.meta.lcon
    #@show x, nlp.meta.lvar, nlp.meta.uvar, mult_x_L, mult_x_U, mult_g

    # doesn't allow inequality constraints.
    comp_vec = [mult_x_L .* max.(0.0,x-nlp.meta.lvar); mult_x_U .* max.(0.0,nlp.meta.uvar-x)]
    comp_vec[[nlp.meta.lvar; nlp.meta.uvar] .== Inf] = 0.0
    comp_vec[[nlp.meta.lvar; nlp.meta.uvar] .== -Inf] = 0.0
    #@show comp_vec
    comp = maximum(comp_vec)
    max_comp = maximum(comp_vec)
    min_comp = minimum(comp_vec)

    fval = obj(nlp,x)

    this_it = OnePhase.generic_alg_history(t,fval,norm_grad_lag,comp,con_vio,y_norm,x_norm)

    push!(hist,this_it)
end

function IPOPT_solver_history(model_builder::Function, solver::IpoptSolver)
    opts = Dict(solver.options);
    tmp_opts = deepcopy(opts)

    temp_m = model_builder()
    nlp = NLPModels.MathProgNLPModel(temp_m)

    hist = Array{OnePhase.generic_alg_history,1}()
    for k = 0:opts[:max_iter]
        m = model_builder()
        tmp_opts[:max_iter] = k
        opts_to_pass = [(val,key) for (val,key) in tmp_opts]
        tmp_solver = IpoptSolver(opts_to_pass)
        setsolver(m, tmp_solver)
        status = solve(m)
        inner = m.internalModel.inner
        add_solver_results!(hist, nlp, inner, k)
        if status != :UserLimit
            break
        end
    end

    return hist
end

function IPOPT_solver_history(nlp::NLPModels.AbstractNLPModel, solver::IpoptSolver; mod=1::Int64) #max_it::Int64; bound_relax_factor::Float64=0.0,tol::Float64=1e-8)
    opts = Dict(solver.options);
    tmp_opts = deepcopy(opts)

    hist = Array{OnePhase.generic_alg_history,1}()
    for k = 0:opts[:max_iter]
        if k % mod == 0
            tmp_opts[:max_iter] = k
            opts_to_pass = [(val,key) for (val,key) in tmp_opts]
            tmp_solver = IpoptSolver(opts_to_pass)
            mp = NLPtoMPB(nlp,tmp_solver)
            status = MathProgBase.optimize!(mp)
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
    for (solver_name,solver) in solver_dic
        if isa(solver,IpoptSolver)
            hist_dic[solver_name] = IPOPT_solver_history(build_nlp, solver);
        else
            m = build_nlp()
            setsolver(m,solver)
            solve(m)
            hist_dic[solver_name] = OnePhase.major_its_only(m.internalModel.inner.hist);
        end
    end
    return hist_dic
end
