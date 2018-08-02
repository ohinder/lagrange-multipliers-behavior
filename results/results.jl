# this file reproduces all the results in the paper
include("../src/shared.jl")

max_it = 100;

solver_dic = Dict(
"One Phase" => OnePhase.OnePhaseSolver(term!max_it= max_it, output_level=0),
# turn off acceptable iterations in IPOPT
"Ipopt" => IpoptSolver(print_level=2, max_iter = max_it, bound_relax_factor=0.0, acceptable_iter=99999, acceptable_tol=Inf, acceptable_compl_inf_tol=Inf, acceptable_constr_viol_tol=Inf, nlp_scaling_method="none"),
"Ipopt w. perturb" => IpoptSolver(print_level=2, max_iter = max_it, acceptable_iter=99999, acceptable_tol=Inf, acceptable_compl_inf_tol=Inf, acceptable_constr_viol_tol=Inf, nlp_scaling_method="none")
)

#####################
## Linear programs ##
#####################
# Figure 1:
# Comparison of the dual variable value (vertical axis) along iterations (horizontal axis) of constraint
# (7d) using IPOPT and a well-behaved IPM [Hinder, 2017] as the perturbation δ is changed.
include("toy_lp.jl")

# Figure 2:
# Comparison on a sample netlib problem
include("sample-netlib.jl")

# Figure 3:
# Comparison of dual multipliers on last 10% of iterations.
include("dual-multipliers.jl")

########################
## Nonconvex programs ##
########################

# Figure 4
# linear complementarity problem
include("comp.jl")

# Figure 5
# Drinking water problem
include("drink.jl")

# Figure 6
# JUNKTURN
include("junk.jl")

######################
## bonus material ####
######################
# TODO
