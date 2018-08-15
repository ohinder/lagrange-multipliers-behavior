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
    sy_inf = OnePhase.get_col(hist,:sy_inf)

    label = ["dual residual" "primal residual" "dual variable" "complementarity"]

    y_dic = Dict(
        "dual residual" => dual_res,
        "primal residual" => primal_res,
        "dual variable" => y_norm,
        "complementarity" => sy_inf
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

function Plot_multiple_solver_trajectories(hist_dic::Dict{String,Array{OnePhase.abstract_alg_history,1}},ylims::Array{Float64,1},display_order::Array{String,1})
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

function Table_from_history(hist_dic::Dict{String,Array{OnePhase.abstract_alg_history,1}},status_dic::Dict{String,Symbol},display_order::Array{String,1};frac::Float64=1.0)
    #df = DataFrame(solver)
    its = Array{Int64,1}()
    max_dual = Array{Float64,1}()
    strict_comp = Array{Float64,1}()
    primal_feas = Array{Float64,1}()
    status_list = Array{Symbol,1}()
    dual_res = Array{Float64,1}()
    primal_res_ls = Array{Float64,1}()
    obj_val_ls = Array{Float64,1}()
    # add sy_inf???

    for solver_name in display_order
        hist = hist_dic[solver_name];

        push!(its,hist[end].t)

        y_norm_hist = OnePhase.get_col(hist,:y_norm);
        push!(max_dual,maximum(slice_frac(y_norm_hist,frac)))

        strict_comp_hist = OnePhase.get_col(hist,:strict_comp);
        push!(strict_comp,minimum(slice_frac(strict_comp_hist,frac)))

        push!(status_list,status_dic[solver_name])

        push!(dual_res,OnePhase.get_col(hist,:norm_grad_lag)[end])

        primal_res = OnePhase.get_col(hist,:primal_residual)[end];
        push!(primal_res_ls,primal_res)

        fval = OnePhase.get_col(hist,:fval)[end];
        push!(obj_val_ls,fval)
    end

    df = DataFrame(
        solvers=display_order,
        iterations=its,
        max_dual=max_dual,
        strict_comp=strict_comp,
        #status_list=status_list,
        dual_res=dual_res,
        primal_res=primal_res_ls,
        fval=obj_val_ls
    )

    return df
end


function table_nl(stream)
    write(stream,"\\\\ \n ")
end

function df_to_latex!(stream,df::DataFrame)
    write(stream,"")

    name_list = names(df)
    for j in 1:length(name_list)
        name = name_list[j]
        spaced_name = replace(String(name),"_"," ")
        write(stream,spaced_name)
        if j != length(name_list)
            write(stream," \& ")
        end
    end
    table_nl(stream)
    write(stream,"\\hline")
    write(stream,"\n")

    for i in 1:size(df, 1)
        for j in 1:size(df,2)
            val = df[i,j]
            if isa(val,Float64) && val != 0.0 && !isinf(val)
                neg = ""
                if val < 0.0
                    val = -val
                    neg = "-"
                end
                pow = floor(Int,log(val)/log(10))
                val = round(val / 10.0^pow,1)
                if pow != 0.0
                    write(stream,"\$ $(neg)$(val) \\times 10\^\{$pow\} \$")
                else
                    write(stream,"\$ $val \$")
                end
            elseif isa(val,Symbol)
                spaced_val = replace(String(val),"_"," ")
                write(stream, spaced_val)
            else
                write(stream,"$val")
            end
            if j != size(df,2)
                write(stream," \& ")
            end
        end
        if i < size(df, 1)
            table_nl(stream)
        end
    end
end

#function df_to_tabular(df::DataFrame,filename::String;caption=""::String)
#    file = open(filename,"w")
#end

function latex_begin_tabular!(stream,df::DataFrame)
    write(stream, "\\begin\{tabular\}\{")
    for j in 1:(size(df,2)-1)
        write(stream,"c ")
    end
    write(stream," c \}")
    write(stream,"\n")
end

function latex_begin_heading!(stream,heading::String,ncol::Int64)
    write(stream, "\\bottomrule\n")
    write(stream,"\\multicolumn\{$ncol\}\{c\}\{$heading\} \\\\ \n")
end

function latex_end_tabular!(stream)
    write(stream, "\\end{tabular} \n")
end

#function add_to_table(stream,df::DataFrame)
#    df_to_latex!(stream,df)
#end
