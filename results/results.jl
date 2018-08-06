################################################################################################
#### This file reproduces all the results in the paper:
#### On the behavior of Lagrange multipliers in convex and non-convex infeasible interior point methods
#### Gabriel Haeser, Oliver Hinder, Yinyu Ye.
#### https://arxiv.org/pdf/1707.07327.pdf
################################################################################################

include("../src/shared.jl")

# check that appropriate folders exist
if !isdir("figures/") mkdir("figures/") end
if !isdir("tables/") mkdir("tables/") end

# keep the same for all solvers
max_it = 100;
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
fig = figure("Trajectory $problem_name",figsize=(9,5))
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

include("circle_example.jl")
hist_dic = Record_solver_histories(solver_dic, circle)
ylims = [1e-6,1e8]
fig = figure("Circle example",figsize=(9,5))
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)
PyPlot.savefig("figures/circle.pdf")
PyPlot.close()

println("************************************************************************************")
println("Computing Figure 5 (complementarity problem) ")
println("************************************************************************************")
println("")

include("comp_example.jl")
hist_dic = Record_solver_histories(solver_dic, simple_comp2)
ylims = [1e-6,1e8]
fig = figure("Complementarity example",figsize=(9,5))
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)

PyPlot.savefig("figures/comp.pdf")
PyPlot.close()

println("************************************************************************************")
println("Computing Figure 6 (drinking water problem)")
println("************************************************************************************")
println("")

include("drink_example.jl")
hist_dic = Record_solver_histories(solver_dic, build_drink5)
ylims = [1e-6,1e8]
fig = figure("Drink example",figsize=(9,5))
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)

PyPlot.savefig("figures/drink.pdf")
PyPlot.close()
