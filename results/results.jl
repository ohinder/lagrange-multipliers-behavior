# this file reproduces all the results in the paper
include("../src/shared.jl")

# keep the same for all solvers
max_it = 3000;
tol = 1e-8
print_level=2

# define solvers
# turn off acceptable iterations in IPOPT
solver_dic = Dict(
"One Phase" => OnePhase.OnePhaseSolver(term!max_it= max_it, output_level=0+print_level),
"Ipopt w/o perturb" => IpoptSolver(tol=tol,print_level=2+print_level, max_iter = max_it, bound_relax_factor=0.0, acceptable_iter=99999, acceptable_tol=tol, acceptable_compl_inf_tol=tol, acceptable_constr_viol_tol=tol, nlp_scaling_method="none"),
"Ipopt w. perturb" => IpoptSolver(tol=tol,print_level=2+print_level, max_iter = max_it, acceptable_iter=99999, acceptable_tol=tol, acceptable_compl_inf_tol=tol, acceptable_constr_viol_tol=tol) #, nlp_scaling_method="none")
)
# order which solvers appear on plots
display_order = ["One Phase", "Ipopt w/o perturb", "Ipopt w. perturb"]

#####################
## Linear programs ##
#####################
# Figure 1:
# Comparison of the dual variable value (vertical axis) along iterations (horizontal axis) of constraint
# (7d) using IPOPT and a well-behaved IPM [Hinder, 2017] as the perturbation δ is changed.
include("toy_lp.jl")

# Figure 2:
# Comparison on a sample netlib problem
problem_name = "ADLITTLE" # which problem???
LP = read_lp(problem_name,"../netlib"); # DOWNLOAD netlib LP problem
build_LP() = build_LP_model_as_NLP(LP)

hist_dic = Record_solver_histories(solver_dic, build_LP)
ylims = [10.0^(-8.0),10.0^(11.0)]
fig = figure("Trajectory $problem_name",figsize=(9,5))
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)
PyPlot.savefig("figures/trajectory_$problem_name.pdf")
PyPlot.close()

# Figure 3:
# Comparison of dual multipliers on last 10% of iterations.
include("dual-multipliers.jl")

########################
## Nonconvex programs ##
########################

# Figure 4
# Drinking water problem
include("circle_example.jl")
hist_dic = Record_solver_histories(solver_dic, circle)
ylims = [1e-6,1e8]
fig = figure("Circle example",figsize=(9,5))
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)
PyPlot.savefig("figures/circle.pdf")
PyPlot.close()

# Figure 5
# linear complementarity problem
include("comp_example.jl")
hist_dic = Record_solver_histories(solver_dic, simple_comp2)
ylims = [1e-6,1e8]
fig = figure("Complementarity example",figsize=(9,5))
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)

PyPlot.savefig("figures/comp.pdf")
PyPlot.close()

# Figure 6
# Drinking water problem
include("drink_example.jl")
hist_dic = Record_solver_histories(solver_dic, build_drink5)
ylims = [1e-6,1e8]
fig = figure("Drink example",figsize=(9,5))
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)

PyPlot.savefig("figures/drink.pdf")
PyPlot.close()

# Figure 7
# JUNKTURN
include("disks.jl")



######################
## bonus material ####
######################
# TODO
