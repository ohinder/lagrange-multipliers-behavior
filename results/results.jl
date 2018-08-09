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
max_it = 3000;
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
    # turn off `acceptable' termination criteron in IPOPT
    :acceptable_iter=>99999,:acceptable_tol=>tol,:acceptable_compl_inf_tol=>tol,:acceptable_constr_viol_tol=>tol
)
options_Ipopt_without_perturb = deepcopy(options_Ipopt_with_perturb)
options_Ipopt_without_perturb[:bound_relax_factor] = 0.0 # this turns off perturbations of the constraints

# define solvers and options in a dictionary
solver_dic = Dict(
"One Phase" => Dict(
    "solver" => :OnePhase,
    "options" => Dict(:term!max_it => max_it, :output_level=>output_level),
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
#####################
## Linear programs ##
#####################

println("************************************************************************************")
println("Computing Figure 1 (toy linear program) ")
println("************************************************************************************")
println("")
include("toy_lp.jl")

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
println("Computing Figure 2 (sample netlib problem) ")
println("************************************************************************************")
println("")
problem_name = "ADLITTLE" # which problem???
LP = read_lp(problem_name,"../netlib"); # DOWNLOAD netlib LP problem
build_LP() = build_LP_model_as_NLP(LP)

hist_dic = Record_solver_histories(solver_dic, build_LP)
ylims = [10.0^(-8.0),10.0^(11.0)]
fig = figure("Trajectory $problem_name",figsize=figsize)
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)
PyPlot.savefig("figures/trajectory_$problem_name.pdf")
PyPlot.close()

########################
## Nonconvex problems ##
########################

println("************************************************************************************")
println("Computing Figure 4 (circle example) ")
println("************************************************************************************")
println("")

function circle()
    m = Model()
    @variable(m, x[1:2])

    @NLobjective(m, Min, x[1])
    @NLconstraint(m, x[1]^2 + x[2]^2 == 1.0)
    @NLconstraint(m, (x[1]-2.0)^2 + x[2]^2 == 1.0)

    return m
end

hist_dic = Record_solver_histories(solver_dic, circle)
ylims = [1e-6,1e8]
fig = figure("Circle example",figsize=figsize)
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)
PyPlot.savefig("figures/circle.pdf")
PyPlot.close()

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

hist_dic = Record_solver_histories(solver_dic, simple_comp)
ylims = [1e-6,1e8]
fig = figure("Complementarity example",figsize=figsize)
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)

PyPlot.savefig("figures/comp.pdf")
PyPlot.close()

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

hist_dic = Record_solver_histories(solver_dic, build_drink)
ylims = [1e-6,1e10]
fig = figure("Drink example",figsize=figsize)
#add_duals_to!(df_duals,"Drink",hist_dic,:y_norm)
#add_strict_comp_to!(df_comp,"Drink",hist_dic,:strict_comp)
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)
#write.csv(df_duals)
PyPlot.savefig("figures/drink.pdf")
PyPlot.close()
