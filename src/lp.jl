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
    setsolver(mod, GurobiSolver(FeasibilityTol=1e-9,BarHomogeneous=1))
    JuMP.build(mod)
    status = solve(mod)
    string_status = (status == :Optimal) ? "true" : "false"
    return string_status
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

function test_solver_on_LP1(LP::LP_problem_data, solver; return_sol=true)
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
end
