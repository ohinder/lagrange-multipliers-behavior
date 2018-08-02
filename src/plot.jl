using Plots, PyPlot

function hist_plot_duals(solver_duals::Dict{String,Dict{String,Float64}})
    duals_IP = collect(values(solver_duals["Ipopt"]));
    duals_OP = collect(values(solver_duals["One Phase"]));
    bins = 0:12
    histogram([log.(duals_IP)/log(10), log.(duals_OP)/log(10)],bins=bins,fillalpha=0.8,xlabel="log_10(max_dual)",label=["Ipopt", "One Phase"],ylabel="number")
end

function plot_duals(solver_duals::Dict{String,Dict{String,Float64}})
    duals_IP = sort(collect(values(solver_duals["Ipopt"])));
    duals_OP = sort(collect(values(solver_duals["One Phase"])));
    prop = collect(1:length(duals_OP))/length(duals_OP)
    plot(prop,[duals_IP, duals_OP],xlabel="proportion of problems",ylabel="maximum dual variable value",label=["Ipopt", "One Phase"],yscale=:log10,legend=:left)
end

function PyPlot_duals(solver_duals::Dict{String,Dict{String,Float64}},ls::Dict{String,String})
    #duals_IP = sort(collect(values(solver_duals["Ipopt"])));
    #duals_OP = sort(collect(values(solver_duals["One Phase"])));

    for (label, data_dic) in solver_duals
        data = sort(collect(values(data_dic)))
        prop = collect(1:length(data))/length(data)
        semilogy(prop, data, label=label, color="black", linestyle=ls[label])
    end
    xlabel("proportion of problems")
    #, xlabel="proportion of problems",ylabel="maximum dual variable value",label=["Ipopt", "One Phase"],yscale=:log10,legend=:left)
end

#=function hist_plot_duals_different_plots(solver_duals::Dict{String,Dict{String,Float64}})
    duals_IP = log.(collect(values(solver_duals["Ipopt"])))/log(10);
    duals_OP = log.(collect(values(solver_duals["OnePhase"])))/log(10);
    combined = [duals_OP; duals_IP]
    bins = 0:12
    ylim = (minimum(combined),maximum(combined))
    histogram([duals_IP,duals_OP],layout=2,bins=bins,xlabel="log_10(max_dual)",legend = false,title=["Ipopt" "One Phase"],ylabel="number of problems",ylim=ylim)
    #histogram(duals_OP,bins=bins,xlabel="log_10(max_dual)",label="One Phase",ylabel="number",ylim=ylim)
end=#

function Plot_jl_trajectory(x, y, label, ylim)
    if ylim == false
        min_y = minimum(y[:])
        max_y = maximum(y[:])
    else
        min_y = ylim[1]
        max_y = ylim[2]
    end
    min_log_y = round(log(min_y)/log(10)+1)
    max_log_y = round(log(max_y)/log(10)-1)
    xticks = 1:10:length(y[:,1])
    yticks = (10.0).^(min_log_y:2:max_log_y)

    Plots.plot(x, y, label=label, legend=:bottom, xticks=xticks, xlabel="# it", yscale = :log10, yticks=yticks, ylim=(min_y,max_y), ylabel="valve")
end

function plot_trajectory(hist; ylim=false,pkg=:Plot_jl)
    assert(length(hist) > 0)

    dual = OnePhase.get_col(hist,:norm_grad_lag)
    primal = OnePhase.get_col(hist,:con_vio)
    y_norm = OnePhase.get_col(hist,:y_norm)
    comp = OnePhase.get_col(hist,:comp)

    #scaled_dual = dual ./ (y_norm + 1.0)

    label = ["dual" "primal" "max dual" "comp"]
    if pkg == :Plot_ji
        y = [dual primal y_norm comp]; # These are the plotting data
        x = 0:(length(y[:,1])-1);
        Plot_jl_trajectory(x, y, label, ylim)
    elseif pkg == :PyPlot
        y_dic = Dict(
            "dual" => dual,
            "primal" => primal,
            "max dual" => y_norm,
            "comp" => comp
        )

        line_style = Dict(
            "dual" => "-",
            "primal" => "-.",
            "max dual" => ":",
            "comp" => "--"
        )
        x = 0:(length(dual)-1)
        PyPlot_trajectory(x, y_dic, line_style)

        if ylim != false
            ax = gca()
            ax[:set_ylim](ylim)
        end
    else
        error("")
    end
end

function PyPlot_trajectory(x, y_dic::Dict{String,Array{Float64,1}}, ls::Dict{String,String})
    for (label, data) in y_dic
        semilogy(x, data, color="black", linestyle=ls[label], label="$label", basey=10)
    end
end
