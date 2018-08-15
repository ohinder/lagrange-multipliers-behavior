# load junk and plot dual multipliers of IPOPT and one phase
# also point out that the final point is relatively close.


using CUTEst

#
#nlp = CUTEstModel("HVYCRASH")
#nlp = CUTEstModel("JUNKTURN") # huge duals
nlp = CUTEstModel("DISCS")



cutest_IPOPT_hist, cutest_status_dic = IPOPT_solver_history(nlp, tmp_solver)

my_pars = OnePhase.create_pars_JuMP(solver_dic["One Phase"].options)
iter, status, full_OP_hist, t, err, timer = OnePhase.one_phase_solve(nlp,my_pars);
OP_hist = OnePhase.major_its_only(full_OP_hist);

#x = OnePhase.get_original_x(iter)

ylims = [1e-7,1e10]
Plot_multiple_solver_dual_histories(hist_dic,ylims,display_order)

PyPlot.savefig("figures/disks.pdf")
PyPlot.close()
