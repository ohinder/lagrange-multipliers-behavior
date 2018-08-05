####################################
### PLOT A TOY PROBLEM
####################################
using JuMP, Ipopt, PyPlot

########################
## code for toy model ##
########################
function build_toy1(perturb::Float64)
    # from previous paper
    m = Model()
    @variable(m, x)
    @variable(m, s[1:2] >= -perturb)

    @NLobjective(m, Min, 0.0 )
    @constraint(m, x + s[2] == 1.0 )
    @constraint(m, x - s[1] == 1.0 )

    return m
end

function build_toy2(perturb::Float64)
    # new paper,
    m = Model()
    @variable(m, x)

    @NLobjective(m, Min, 0.0 )
    @constraint(m, x <= 1.0 + perturb )
    @constraint(m, x >= 1.0 - perturb )

    return m
end

build_toy = build_toy2; # which model should we use.

ls = Dict(1e-3 => "-", 1e-6 => "-.", 1e-9 => "--", 0.0 => ":") # line styles for plot
ylims = [1e-6, 1e10] # y-axis maximum and min value
delta_set = [1e-3, 1e-6, 1e-9, 0.0] # choice of perturbations

###########
## IPOPT ##
###########
subplot(121)
for delta = delta_set
  build_LP() = build_toy(delta)
  hist = IPOPT_solver_history(build_LP, solver_dic["Ipopt"])
  dual_hist = OnePhase.get_col(hist, :y_norm)

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
for delta = delta_set
  m = build_toy(delta)
  setsolver(m,solver_dic["One Phase"])
  status = solve(m)
  dual_hist = OnePhase.get_col(m.internalModel.inner.hist, :y_norm)
  semilogy(1:length(dual_hist), dual_hist, color="black", linestyle=ls[delta], label="δ = $delta", basey=10)
end

title("Well behaved IPM")
xlabel("iteration")
ax = gca()
ax[:set_ylim](ylims)

PyPlot.savefig("figures/toy_lp.pdf")

PyPlot.close()
