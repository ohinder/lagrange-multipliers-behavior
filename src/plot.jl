using PyPlot

function Plot_duals(solver_duals::Dict{String,Dict{String,Float64}},ls::Dict{String,String}; ylim=false)
    for (label, data_dic) in solver_duals
        data = sort(collect(values(data_dic)))
        prop = collect(1:length(data))/length(data)
        semilogy(prop, data, label=label, color="black", linestyle=ls[label])
    end
    xlabel("proportion of problems")

    if ylim != false
        ax = gca()
        ax[:set_ylim](ylim)
    end
end

function Plot_trajectory(hist; ylim=false)
    assert(length(hist) > 0)

    its = OnePhase.get_col(hist,:t)
    dual_res = OnePhase.get_col(hist,:norm_grad_lag)
    primal_res = OnePhase.get_col(hist,:primal_residual)
    y_norm = OnePhase.get_col(hist,:y_norm)
    comp = OnePhase.get_col(hist,:comp)

    label = ["dual residual" "primal residual" "dual variable" "complementarity"]

    y_dic = Dict(
        "dual residual" => dual_res,
        "primal residual" => primal_res,
        "dual variable" => y_norm,
        "complementarity" => comp
    )

    line_style = Dict(
        "dual residual" => "-.",
        "primal residual" => "--",
        "dual variable" => "-",
        "complementarity" => ":"
    )
    Plot_trajectory(its, y_dic, line_style)

    if ylim != false
        ax = gca()
        ax[:set_ylim](ylim)
    end

    xlabel("iteration")
end

function Plot_trajectory(x, y_dic::Dict{String,Array{Float64,1}}, ls::Dict{String,String})
    for (label, data) in y_dic
        semilogy(x, data, color="black", linestyle=ls[label], label="$label", basey=10)
    end
end

function Plot_multiple_solver_dual_histories(hist_dic::Dict{String,Array{OnePhase.abstract_alg_history,1}},ylims::Array{Float64,1},display_order::Array{String,1})
    @assert(length(ylims) == 2)
    n = length(hist_dic)
    i = 0
    for solver_name in display_order
        hist = hist_dic[solver_name]
        i += 1
        PyPlot.subplot(1,n,i)
        Plot_trajectory(hist,ylim=ylims)
        title(solver_name)

        if i == 1
            ylabel("L-infinity norm")
            legend()
        else
            ax = gca()
            ax[:set_yticklabels]([])
        end
    end
end
