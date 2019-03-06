################################################################################################
#### This file reproduces all the results in the paper:
#### On the behavior of Lagrange multipliers in convex and non-convex infeasible interior point methods
#### Gabriel Haeser, Oliver Hinder, Yinyu Ye.
#### https://arxiv.org/pdf/1707.07327.pdf
################################################################################################

include("../src/shared.jl") # functions that we use in results.jl

# check that appropriate folders exist
if !isdir("figures/") mkdir("figures/") end
if !isdir("tables/") mkdir("tables/") end

# keep the same for all solvers
max_it = 300;
tol = 1e-6;
output_level=0

println("Solver parameters (shared by all solvers):")
@show max_it, tol, output_level

# define options for IPOPT
# see https://www.coin-or.org/Ipopt/documentation/node40.html
options_Ipopt_with_perturb = Dict(
    :tol=>tol, # termination tolerance
    :print_level=>2+output_level,
    :max_iter => max_it,
    #:nlp_scaling_method => "none",
    # turn off `acceptable' termination criteron in IPOPT
    :acceptable_iter=>99999,:acceptable_tol=>tol,:acceptable_compl_inf_tol=>tol,:acceptable_constr_viol_tol=>tol
)

options_OnePhase = Dict(
    :term!max_it => max_it,
    :output_level => output_level,
    :term!tol_opt => tol
);
options_Ipopt_without_perturb = deepcopy(options_Ipopt_with_perturb)
options_Ipopt_without_perturb[:bound_relax_factor] = 0.0 # this turns off perturbations of the constraints

# define solvers and options in a dictionary
solver_dic = Dict(
"One Phase" => Dict(
    "solver" => :OnePhase,
    "options" => options_OnePhase,
),
"Ipopt w/o perturb" => Dict(
    "solver" => :Ipopt,
    "options" => options_Ipopt_without_perturb,
),
"Ipopt w. perturb" => Dict(
    "solver" => :Ipopt,
    "options" => options_Ipopt_with_perturb,
));

# order which solvers appear on plots
display_order = ["One Phase", "Ipopt w/o perturb", "Ipopt w. perturb"]

figsize = (9,4);

###############################
## write the options to a file
###############################
function opts_to_df(dic::Dict) return DataFrame(option=collect(keys(dic)),value=collect(values(dic))) end

df_options_OnePhase = opts_to_df(options_OnePhase)
df_options_Ipopt_without_perturb = opts_to_df(options_Ipopt_without_perturb)
df_options_Ipopt_with_perturb = opts_to_df(options_Ipopt_with_perturb)

latex_tbl = open("tables/options_table.txt","w")
latex_begin_tabular!(latex_tbl,df_options_OnePhase)
latex_begin_heading!(latex_tbl,"Ipopt",2)
df_to_latex!(latex_tbl, df_options_Ipopt_without_perturb)
write(latex_tbl,"\\\\ \n")
latex_begin_heading!(latex_tbl,"One Phase",2)
df_to_latex!(latex_tbl, df_options_OnePhase)
write(latex_tbl,"\\\\ \n")
latex_end_tabular!(latex_tbl)
close(latex_tbl)

#####################
## Linear programs ##
#####################

println("************************************************************************************")
println("Computing Figure 3 (Comparision of dual multipliers on netlib)")
println("************************************************************************************")
println("")
include("dual-multipliers.jl")

# For the remainder of the figures lets have a look at the output
solver_dic["One Phase"]["options"][:output_level] = 2;
solver_dic["Ipopt w/o perturb"]["options"][:print_level] = 5;
solver_dic["Ipopt w. perturb"]["options"][:print_level] = 5;

println("************************************************************************************")
println("Computing Figure 1 (toy linear program) ")
println("************************************************************************************")
println("")
include("toy_lp.jl")

println("************************************************************************************")
println("Computing Figure 2 (sample netlib problem) ")
println("************************************************************************************")
println("")
#problem_name = "ADLITTLE" # which problem???
#problem_name = "ISRAEL" # which problem???
problem_name = "BANDM" # which problem???

LP = read_lp(problem_name,"../netlib"); # DOWNLOAD netlib LP problem
build_LP() = build_LP_model_as_NLP(LP)

example_hist_dic, example_status_dic = Record_solver_histories(solver_dic, build_LP, build_LP);
ylims = [10.0^(-8.0),10.0^(12.0)]
df_netlib_example = Table_from_history(example_hist_dic,example_status_dic,display_order,frac=0.2)
CSV.write("tables/netlib_example_table.csv",df_netlib_example)
fig = figure("Trajectory $problem_name",figsize=figsize)
Plot_multiple_solver_trajectories(example_hist_dic,ylims,display_order)

for i = 1:3
    println("Duals for ",display_order[i], " are")
    println(OnePhase.get_col(example_hist_dic[display_order[i]],:y_norm))
end

PyPlot.savefig("figures/trajectory_$problem_name.pdf")
PyPlot.close()

########################
## Nonconvex problems ##
########################

println("************************************************************************************")
println("Computing Figure 4 (circle example) ")
println("************************************************************************************")
println("")

function circle_ineq()
    m = Model()
    @variable(m, x[1:2])

    @NLobjective(m, Min, -(x[1]-1.0)^2 + x[2]^2)
    @NLconstraint(m, x[1]^2 + x[2]^2 <= 1.0)
    @NLconstraint(m, (x[1]-2.0)^2 + x[2]^2 <= 1.0)

    return m
end

function circle_ineq_ipopt()
    m = Model()
    @variable(m, x[1:2])
    @variable(m, s[1:2] >= 0.0)

    @NLobjective(m, Min, -(x[1]-1.0)^2 + x[2]^2)
    @NLconstraint(m, x[1]^2 + x[2]^2 + s[1] == 1.0)
    @NLconstraint(m, (x[1]-2.0)^2 + x[2]^2 + s[2] == 1.0)

    return m
end

circle_hist_dic,circle_status_dic = Record_solver_histories(solver_dic, circle_ineq, circle_ineq_ipopt);
#drink_status_dic["Ipopt w. perturb"] = :eq_mult_error # for this problem Ipopt incorrectly declares optimality -- see console output
#drink_status_dic["Ipopt w/o perturb"] = :eq_mult_error # for this problem Ipopt incorrectly declares optimality -- see console output

ylims = [1e-6,1e11]
fig = figure("Circle example",figsize=figsize)
Plot_multiple_solver_trajectories(circle_hist_dic,ylims,display_order)
PyPlot.savefig("figures/circle.pdf")
PyPlot.close()

df_circle = Table_from_history(circle_hist_dic,circle_status_dic,display_order,frac=0.2)
CSV.write("tables/circle_table.csv",df_circle)

# latex table (open file)
latex_tbl = open("tables/latex_table.txt","w")
latex_ncol = length(names(df_circle))
latex_begin_tabular!(latex_tbl,df_circle)
latex_begin_heading!(latex_tbl,"Intersection of two circles",latex_ncol)
df_to_latex!(latex_tbl, df_circle)
write(latex_tbl,"\\\\ \n")

println("************************************************************************************")
println("Computing Figure 5 (complementarity problem) ")
println("************************************************************************************")
println("")

function simple_comp()
    m = Model()
    @variable(m, x[1:2] >= 0.0)
    #@variable(m, 0.5 >= y >= 0.0)

    @objective(m, Min, 3 * x[1] + x[2])

    @constraint(m, x[1] - 3.0 * x[2] >= 2.0)
    @NLconstraint(m, x[1] * x[2] <= 0.0)

    return m
end

function simple_comp_ipopt()
    m = Model()
    @variable(m, x[1:2] >= 0.0)
    @variable(m, s[1:2] >= 0.0)

    @objective(m, Min, 3 * x[1] + x[2])

    @constraint(m, x[1] - 3.0 * x[2] - s[1] == 2.0)
    @NLconstraint(m, x[1] * x[2] + s[2] == 0.0)

    return m
end

comp_hist_dic, comp_status_dic = Record_solver_histories(solver_dic, simple_comp, simple_comp_ipopt);
ylims = [1e-6,1e11]
fig = figure("Complementarity example",figsize=figsize)
Plot_multiple_solver_trajectories(comp_hist_dic,ylims,display_order)
PyPlot.savefig("figures/comp.pdf")
PyPlot.close()

df_comp = Table_from_history(comp_hist_dic,comp_status_dic,display_order,frac=0.2)
CSV.write("tables/comp_table.csv",df_comp)

# latex table
latex_begin_heading!(latex_tbl,"Linear program with complementarity constraints",latex_ncol)
df_to_latex!(latex_tbl,df_comp)
write(latex_tbl,"\\\\ \n")



println("************************************************************************************")
println("Computing Figure 6 (drinking water problem)")
println("************************************************************************************")
println("")
##
## on this problem Ipopt fails because of the error
## `Cannot recompute multipliers for feasibility problem.  Error in eq_mult_calculator`

function build_drink()
    m = Model()
    #@variable(m, x[1:3] >= 0.0, start=1.0 )
    #@variable(m, y[1:3], start=1.0)
    #@variable(m, z >= 0.0, start=100.0)
    @variable(m, x[1:3] >= 0.0) #, start=1.0)
    @variable(m, h[1:3] >= 0.0) #, start=1.0)

    @NLobjective(m, Min, h[1] )
    # flow equations
    @constraint(m, x[1] + x[2] == 1.0) # node 1 supplies one unit of water
    @constraint(m, x[1] + x[3] == 0.5) # demand at node 1 must be met
    @constraint(m, x[2] == 0.5) # demand at node 2 must be met

    theta = 1.8

    # pressure equations
    @NLconstraint(m, x[1]^theta - h[1] +  h[2] == 0.0) # pressure loss from node 1 to 2
    @NLconstraint(m, x[2]^theta - h[1] + h[3] == 0.0)  # pressure loss from node 1 to 3
    @NLconstraint(m, x[3]^theta - h[2] + h[3] == 0.0) # pressure loss from node 2 to 3

    return m
end

drink_hist_dic, drink_status_dic = Record_solver_histories(solver_dic, build_drink, build_drink);
drink_status_dic["Ipopt w. perturb"] = :eq_mult_error # for this problem Ipopt incorrectly declares optimality -- see console output
drink_status_dic["Ipopt w/o perturb"] = :eq_mult_error # for this problem Ipopt incorrectly declares optimality -- see console output

ylims = [1e-6,1e11]
fig = figure("Drink example",figsize=figsize)
#add_duals_to!(df_duals,"Drink",hist_dic,:y_norm)
#add_strict_comp_to!(df_comp,"Drink",hist_dic,:strict_comp)
Plot_multiple_solver_trajectories(drink_hist_dic,ylims,display_order)
PyPlot.savefig("figures/drink.pdf")
PyPlot.close()

df_drink = Table_from_history(drink_hist_dic,drink_status_dic,display_order,frac=0.2)
CSV.write("tables/drink_table.csv",df_drink)

# latex table
latex_begin_heading!(latex_tbl,"Drinking water network optimization",latex_ncol)
df_to_latex!(latex_tbl,df_drink)
write(latex_tbl,"\\\\ \n \\bottomrule\n")
latex_end_tabular!(latex_tbl)
close(latex_tbl)
