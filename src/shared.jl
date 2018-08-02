using JuMP, MAT, Ipopt, MathProgBase, NLPModels, CUTEst
using Gurobi # to check if LP is feasible or not
#using OnePhase
include("../../one-phase-2.0/src/OnePhase.jl")
include("lp.jl")
include("plot.jl")

function add_solver_results!(hist::Array{OnePhase.generic_alg_history,1}, nlp::AbstractNLPModel, inner)
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

    t = length(hist)

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
        add_solver_results!(hist, nlp, inner)
        if status == :Optimal
            break
        end
    end

    return hist
end

function IPOPT_solver_history(nlp::CUTEstModel, max_it::Int64; bound_relax_factor::Float64=0.0,tol::Float64=1e-8)
    optimal_it = max_it
    for j = 1:max_it
        mp = NLPtoMPB(nlp, IPOPT_k_its_normal(j,bound_relax_factor,tol));
        status = MathProgBase.optimize!(mp)
        @show status

        if status != -1
            optimal_it = j
            break
        end
    end

    @show optimal_it

    hist = Array{OnePhase.generic_alg_history,1}()
    for k = 1:optimal_it
        mp = NLPtoMPB(nlp, IPOPT_k_its_normal(k,bound_relax_factor,tol));
        status = MathProgBase.optimize!(mp)

        add_solver_results!(hist, nlp, mp.inner)
    end

    return hist
end
