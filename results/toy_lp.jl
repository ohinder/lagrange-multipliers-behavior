####################################
### PLOT A TOY PROBLEM
####################################

########################
## code for toy model ##
########################
function build_toy_ipopt(perturb::Float64)
    m = Model()
    @variable(m, x)
    @variable(m, s[1:2] >= 0.0) # add slacks manually

    @NLobjective(m, Min, 0.0 )
    @constraint(m, x + s[1] == 1.0 + perturb )
    @constraint(m, x - s[2] == 1.0 - perturb )

    return m
end

function build_toy(perturb::Float64)
    # new paper,
    m = Model()
    @variable(m, x)

    @NLobjective(m, Min, 0.0 )
    @constraint(m, x <= 1.0 + perturb )
    @constraint(m, x >= 1.0 - perturb )

    return m
end

#=function confused(perturb::Float64)
    # new paper,
    m = Model()
    @variable(m, x)

    @NLobjective(m, Min, x^2 )
    @constraint(m, x <= 1.0 + perturb )
    @constraint(m, x >= 1.0 - perturb )

    return m
end=#

#build_toy = build_toy2; # which model should we use.

ls = Dict(1e-2 => ":", 1e-5 => "--", 1e-8 => "-.", 0.0 => "-") # line styles for plot
ylims = [1e-6, 1e10] # y-axis maximum and min value
delta_set = sort(collect(keys(ls))) # choice of perturbations

###########
## IPOPT ##
###########
subplot(121)
println("IPOPT dual history ...")
for delta = delta_set
  build_LP() = build_toy_ipopt(delta)
  hist, dual_dic = IPOPT_solver_history(build_LP, solver_dic["Ipopt w/o perturb"])
  dual_hist = OnePhase.get_col(hist, :y_norm)
  @show delta
  @show dual_hist

  semilogy(1:length(dual_hist), dual_hist, color="black", linestyle=ls[delta], label="δ = $delta", basey=10)
end

title("IPOPT")
ax = gca()
ax[:set_ylim](ylims)
xlabel("iteration")
ylabel("maximum dual variable")

legend()

###############
## one-phase ##
###############
subplot(122)
println("OnePhase dual history ...")
for delta = delta_set
  m = build_toy(delta)
  #m = confused(delta)
  setsolver(m,build_solver(solver_dic["One Phase"]))
  status = solve(m)
  hist = OnePhase.major_its_only(m.internalModel.inner.hist)
  dual_hist = OnePhase.get_col(hist, :y_norm)
  @show delta
  @show dual_hist

  semilogy(1:length(dual_hist), dual_hist, color="black", linestyle=ls[delta], label="δ = $delta", basey=10)
end

title("Well behaved IPM")
xlabel("iteration")
ax = gca()
ax[:set_ylim](ylims)

PyPlot.savefig("figures/toy_lp.pdf")

PyPlot.close()
