mutable struct LP_problem_data
    A
    b::Array{Float64,1}
    c::Array{Float64,1}
    lbounds::Array{Float64,1}
    ubounds::Array{Float64,1}
end

mutable struct LP_solution
    x::Array{Float64,1}
    y::Array{Float64,1}
    s::Array{Float64,1}
end

function read_lp(name::String, dir::String)
  lp_data = matopen("$(dir)/$(name).mat")

  A = read(lp_data,"A")
  b = read(lp_data,"b")[:]
  c = read(lp_data,"c")[:]
  lbounds = nothing
  ubounds = nothing

  try
    lbounds = read(lp_data,"lbounds")[:]
    ubounds = read(lp_data,"ubounds")[:]
  catch(e)
    lbounds = zeros(length(c))
    ubounds = Inf * ones(length(c))
  end

  ubounds[ubounds .> 10^10.0] = Inf
  lbounds[lbounds .< -10^10.0] = -Inf

  close(lp_data)

  return LP_problem_data(A, b, c, lbounds, ubounds)
end

function load_netlib(num_nz::Int64, dir::String)
    file_list = readdir("../netlib")
    problem_list = Array{String,1}()
    for file_name in file_list
        if file_name[(end-3):end] == ".mat" #&& filesize("netlib/$file_name") < file_size
          problem_name = file_name[1:(end-4)]
          println("loading...", problem_name)
          LP = read_lp(problem_name, dir)
          if nnz(LP.A) <= num_nz
            push!(problem_list, problem_name)
          end
        end
    end

    return problem_list
end

function build_LP_model(LP::LP_problem_data)
    m = Model()
    @variable(m, LP.lbounds[i] <= x[i=1:length(LP.c)] <= LP.ubounds[i], start = 0.0)
    @objective(m, Min, LP.c' * x)
    equality_con = @constraint(m, LP.A * x .== LP.b )
    return m;
end

function build_LP_model_as_NLP(LP::LP_problem_data)
    m = Model()
    @variable(m, LP.lbounds[i] <= x[i=1:length(LP.c)] <= LP.ubounds[i], start = 0.0)
    @NLobjective(m, Min, sum(LP.c[i] * x[i] for i=1:length(LP.c)))
    equality_con = @constraint(m, LP.A * x .== LP.b )
    return m;
end

#=function test_solver_on_LP1(LP::LP_problem_data, solver; return_sol=true)
  m = Model(solver = solver)

  @variable(m, LP.lbounds[i] <= x[i=1:length(LP.c)] <= LP.ubounds[i])

  @objective(m, Min, LP.c' * x)

  equality_con = @constraint(m, LP.A * x .== LP.b )

  status = JuMP.solve(m)

  #@show status
  if return_sol
      x_val = getvalue(x)

      y_val = getdual(equality_con)
      rc_val = getdual(x)

      return LP_solution(x_val, y_val, rc_val), status
    else
        status
    end
end=#

function perturb_LP!(LP::LP_problem_data,perturb::Float64)
    LP.lbounds -= perturb
    LP.ubounds += perturb
end

function lp_feasible(LP::LP_problem_data, tol::Float64)
    # check if is an LP is feasible after perturbation
    @assert(tol >= 1e-8)
    LP_tmp = deepcopy(LP)
    perturb_LP!(LP_tmp,-tol)
    mod = build_LP_model(LP_tmp)
    setsolver(mod, GurobiSolver(FeasibilityTol=1e-9,BarHomogeneous=1,OutputFlag=0))
    JuMP.build(mod)
    status = solve(mod)
    return status == :Optimal
end

function problems_with_no_interior(problem_list::Array{String,1}, tol_set::Array{Float64,1})
    tol_set = [1e-4, 1e-6, 1e-8]
    status_list = Dict{Float64,Array{Bool,1}}()

    for tol in tol_set
        status_list[tol] = Array{Float64,1}()
        for problem_name in problem_list
          LP = read_lp(problem_name,"../netlib/")
          status = lp_feasible(LP, tol)
          push!(status_list[tol], status)
        end
    end

    return status_list
end

function compare_across_tols(int_status_list::Dict{Float64,Array{Bool,1}},tol_set::Array{Float64,1})
    interior = Dict()
    no_interior = Dict()

    println("tol | # int | # no int")
    for tol in tol_set
      interior[tol] = sum(int_status_list[tol])
      no_interior[tol] = sum(.!int_status_list[tol])
      println(tol, " | ", interior[tol], " | ", no_interior[tol])
    end

    return interior, no_interior
end


function print_solvers_failures(solver_status_dic::Dict{String,Dict{String,Symbol}})
    println("===== Which problems do the solvers fail on? ======")

    for (solver_name, solver_status) in solver_status_dic
        num_failures = 0
        println("==== Solver: ", solver_name, " ======")
        for (prob_name,status) in solver_status
            if status != :Optimal
                println(prob_name, " ", status)
                num_failures += 1;
            end
        end
        println("total failures = ",num_failures)
    end
end

function compute_solver_status(solver_dic::Dict,problem_list::Array{String,1})
    solver_status_dic = Dict{String,Dict{String,Symbol}}()
    for (solver_name,solver) in solver_dic
        solver_status_dic[solver_name] = Dict{String,Symbol}()

        for problem_name in problem_list
          println("Solving $problem_name")
          LP = read_lp(problem_name,"../netlib/")
          m = build_LP_model_as_NLP(LP)
          setsolver(m, solver)
          status = solve(m)
          solver_status_dic[solver_name][problem_name] = status
        end
    end

    return solver_status_dic
end

function all_solvers_succeed(solver_status_dic::Dict{String,Dict{String,Symbol}},problem_list::Array{String,1})
    # construct a list of problems for which all solvers succeed
    both_succeed = Array{String,1}();
    for prob_name in problem_list
        succeed = true
        for (solver_name,solver_status) in solver_status_dic
            if solver_status[prob_name] != :Optimal
                succeed = false
            end
        end
        if succeed == true
            push!(both_succeed,prob_name)
        end
    end

    return both_succeed
end

function get_mult_inf_norm(inn::OnePhase.OnePhaseProblem)
    return norm(inn.lambda,Inf)
end

function get_mult_inf_norm(inn::IpoptProblem)
    return max(norm(inn.mult_g,Inf),norm(inn.mult_x_L,Inf),norm(inn.mult_x_U,Inf))
end

function get_maximum_duals_at_last_iterate(solver_dic::Dict, problem_list::Array{String,1})
    # generate distribution on the **maximum** dual multiplier values at the **last iterate**
    solver_duals = Dict{String,Dict{String,Float64}}()
    for (solver_name,solver) in solver_dic
        solver_duals[solver_name] = Dict{String,Float64}()

        for problem_name in problem_list
          println("Solving $problem_name")
          LP = read_lp(problem_name,"../netlib/")
          m = build_LP_model_as_NLP(LP)
          setsolver(m, solver)
          status = solve(m)
          solver_duals[solver_name][problem_name] = get_mult_inf_norm(m.internalModel.inner)
        end
    end

    return solver_duals
end

function get_maximum_duals_all_iterates(solver_dic::Dict, problem_list::Array{String,1}; frac::Float64=1.0)
    # generate distribution on the **maximum** dual multiplier values at the **last iterate**
    solver_duals = Dict{String,Dict{String,Float64}}()
    for (solver_name,solver) in solver_dic
        solver_duals[solver_name] = Dict{String,Float64}()

        for problem_name in problem_list
          println("Solving $problem_name")
          LP = read_lp(problem_name,"../netlib/")

          if isa(solver, IpoptSolver)
              build_LP() = build_LP_model_as_NLP(LP)
              hist = IPOPT_solver_history(build_LP, solver)
          elseif solver_name == "One Phase"
                m = build_LP_model_as_NLP(LP)
                setsolver(m,OnePhase.OnePhaseSolver(output_level=2))
                solve(m)

                hist = OnePhase.major_its_only(m.internalModel.inner.hist)
          else
              @show solver_name, typeof(solver)
              error("Invalid solver_name")
          end
          vals = OnePhase.get_col(hist,:y_norm)
          if frac == 1.0
              j = 1
          elseif frac >= 0.0 && frac < 1.0
              j = round(Int,ceil(length(vals) * (1.0-frac)))
          else
              error("")
          end
          slice = vals[j:end]
          solver_duals[solver_name][problem_name] = maximum(slice)
        end
    end

    return solver_duals
end
