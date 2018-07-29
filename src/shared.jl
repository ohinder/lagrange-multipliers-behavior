using JuMP, MAT, Ipopt, Plots, MathProgBase, NLPModels, CUTEst
#using OnePhase
include("../../one-phase-2.0/src/OnePhase.jl")
include("lp.jl")

function IPOPT_k_its_acceptable(k::Int64,bound_relax_factor::Float64,tol::Float64)
    return IpoptSolver(tol=tol, print_level=1, acceptable_tol=Inf, acceptable_iter=k, max_iter=3000, acceptable_constr_viol_tol=Inf, acceptable_dual_inf_tol=Inf, acceptable_compl_inf_tol=Inf, acceptable_obj_change_tol=Inf)
end

function IPOPT_k_its_normal(k::Int64,bound_relax_factor::Float64,tol::Float64)
    return IpoptSolver(tol=tol, bound_relax_factor=bound_relax_factor, print_level=1, acceptable_iter=9999999, max_iter=k)
end

function IPOPT_solver_history(model_builder::Function, max_it::Int64; bound_relax_factor::Float64=0.0,tol::Float64=1e-8)
    optimal_it = max_it
    for j = 1:max_it
        m = model_builder()
        setsolver(m, IPOPT_k_its_normal(j,bound_relax_factor,tol))

        status = solve(m)

        if status == :Optimal
            optimal_it = j
            break
        end
    end

    @show optimal_it

    temp_m = model_builder()
    nlp = NLPModels.MathProgNLPModel(temp_m)

    hist = Array{OnePhase.generic_alg_history,1}()
    for k = 1:(optimal_it+1)
        m = model_builder()
        setsolver(m, IPOPT_k_its_normal(k,bound_relax_factor,tol))
        status = solve(m)
        inner = m.internalModel.inner
        OnePhase.add_solver_results!(hist, nlp, inner)
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

        OnePhase.add_solver_results!(hist, nlp, mp.inner)
    end

    return hist
end


function plot_trajectory(x, y, label; min_y=false,max_y=false)
    if min_y == false
        min_y = minimum(y[:])
    end
    if max_y == false
        max_y = maximum(y[:])
    end
    min_log_y = round(log(min_y)/log(10)+1)
    max_log_y = round(log(max_y)/log(10)-1)
    xticks = 1:10:length(y[:,1])
    yticks = (10.0).^(min_log_y:2:max_log_y)

    Plots.plot(x, y, label=label, legend=:bottom, xticks=xticks, xlabel="# it", yscale = :log10, yticks=yticks, ylim=(min_y,max_y), ylabel="valve")
end

function plot_trajectory(hist; min_y=false,max_y=false)
    assert(length(hist) > 0)

    dual = OnePhase.get_col(hist,:norm_grad_lag)
    primal = OnePhase.get_col(hist,:con_vio)
    y_norm = OnePhase.get_col(hist,:y_norm)
    comp = OnePhase.get_col(hist,:comp)

    scaled_dual = dual ./ (y_norm + 1.0)

    y = [dual primal y_norm comp]; # These are the plotting data
    x = 1:length(y[:,1]);

    label = ["dual" "primal" "max dual" "comp"]
    plot_trajectory(x, y, label, min_y=min_y,max_y=max_y)
end
