max_problem_size = 1000;

# calculate which Netlib LP problems have no interior
problem_list = load_netlib(max_problem_size,"../netlib/")
tol_set = [1e-4, 1e-6, 1e-8];
probs_no_int_tol_set = problems_with_no_interior(problem_list, tol_set);

println("Check tolerance doesn't really matter ...")
compare_across_tols(probs_no_int_tol_set, tol_set);
has_int = probs_no_int_tol_set[tol_set[1]];
problems_with_int = problem_list[has_int]; # problems with interior
problems_with_no_int = problem_list[.!has_int]; # problems with no interior
println("=================")

# which problems have an interior? write to CSV.
using DataFrames,CSV
df = DataFrame(problem_name=problem_list, has_int=has_int)
CSV.write("tables/has_int.csv", df)

# compute number of failures for all three methods
solver_status_dic = compute_solver_status(solver_dic, problem_list);
print_solvers_failures(solver_status_dic);

#restricted_problem_list = intersect(problem_list,problems_with_no_int)
problems_where_all_solvers_succeed = all_solvers_succeed(solver_status_dic, problem_list);
all_solvers_succeed_and_no_int = intersect(problems_where_all_solvers_succeed,problems_with_no_int)
all_solvers_succeed_and_int = intersect(problems_where_all_solvers_succeed,problems_with_int)

# generate distribution on the **maximum** dual multiplier values over **all iterates**
solver_duals = get_maximum_duals_all_iterates(solver_dic, problems_where_all_solvers_succeed, frac=1.0)
solver_duals_20_percent = get_maximum_duals_all_iterates(solver_dic, problems_where_all_solvers_succeed,frac=0.2)

line_style = Dict("One Phase" => "--", "Ipopt" => "-", "Ipopt w. perturb" => "-.")

PyPlot_duals(solver_duals,line_style)
legend()
PyPlot.savefig("figures/dual_distribution.pdf")
PyPlot.close()

PyPlot_duals(solver_duals_20_percent,line_style)
legend()
PyPlot.savefig("figures/dual_distribution_20%.pdf")
PyPlot.close()

function restrict_results(res::Dict{String,Dict{String,Float64}}, restricted_problem_list::Array{String,1})
    restricted_res = Dict{String,Dict{String,Float64}}()
    for (solver_name,data) in res
        restricted_res[solver_name] = Dict{String,Float64}()
        for problem_name in restricted_problem_list
            restricted_res[solver_name][problem_name] = data[problem_name]
        end
    end
    return restricted_res
end

# effect is clearer when we look only at problems with no interior
solver_duals_20_percent_no_int = restrict_results(solver_duals_20_percent,all_solvers_succeed_and_no_int)
PyPlot_duals(solver_duals_20_percent_no_int,line_style)
legend()
PyPlot.savefig("figures/dual_distribution_20%_no_int.pdf")
PyPlot.close()

solver_duals_20_percent_int = restrict_results(solver_duals_20_percent,all_solvers_succeed_and_int)
PyPlot_duals(solver_duals_20_percent_int,line_style)
legend()
PyPlot.savefig("figures/dual_distribution_20%_int.pdf")
PyPlot.close()
