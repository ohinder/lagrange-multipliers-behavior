## compute dual variables sequence
problem_name = "ADLITTLE" # which problem???
LP = read_lp(problem_name,"../netlib"); # DOWNLOAD netlib LP problem
build_LP() = build_LP_model_as_NLP(LP)


ylim = (10.0^(-8.0),10.0^(11.0))

## compute trajectory for one-phase
m = build_LP()
setsolver(m,solver_dic["One Phase"])
solve(m)
OP_hist = OnePhase.major_its_only(m.internalModel.inner.hist);

# plot the trajectory for one-phase
subplot(131)
plot_trajectory(OP_hist,ylim=ylim,pkg=:PyPlot)
title("One Phase")
xlabel("iteration")
ylabel("value")

legend() # legend is same for all plots

## because this IPOPT does not provide this information by default we have to use a hack
## which takes significantly longer than just solving the linear program normally
## (we solve the LP with a maximum iteration of k = 1,...,max_iter)
IPOPT_hist = IPOPT_solver_history(build_LP, solver_dic["Ipopt"]);

# PLOT THE TRAJECTORY FOR IPOPT
subplot(132)
plot_trajectory(IPOPT_hist,ylim=ylim,pkg=:PyPlot)
title("IPOPT")
ax = gca()
ax[:set_yticklabels]([])

p_LP = deepcopy(LP)
perturb_LP!(p_LP,1e-8) # ** CHECK THIS MATCHES PAPER!!!! **
p_build_LP() = build_LP_model_as_NLP(p_LP) # add peturbation
IPOPT_hist = IPOPT_solver_history(p_build_LP, solver_dic["Ipopt"]);

# PLOT THE TRAJECTORY FOR IPOPT with perturbation
subplot(133)
plot_trajectory(IPOPT_hist,ylim=ylim,pkg=:PyPlot)
title("IPOPT w. peturb")
ax = gca()
ax[:set_yticklabels]([])

PyPlot.savefig("figures/trajectory_$problem_name.pdf")

PyPlot.close()
