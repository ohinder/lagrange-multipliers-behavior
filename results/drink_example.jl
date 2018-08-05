# drinking water problem?
function build_drink1()
    m = Model()
    @variable(m, 0 <= x <= 2, start = 1 )
    @variable(m, y <= 0)
    @variable(m, z[1:2] >= 0)

    @objective(m, Min, z[1]^2 + z[2]^2 + z[1] * z[2]  )
    @NLconstraint(m, x^2 - y == 0.0)
    @NLconstraint(m, z[1] + z[2]^3 - 1.0 == 0.0)

    return m
end

# drinking water problem?
function build_drink2()
    m = Model()
    @variable(m, x[1:6] >= 0.0,start=1.0 )
    @variable(m, y[1:3])
    @variable(m, z[1:3] >= 0.0)

    @NLobjective(m, Min,
        z[1]^2 + 2.0 * z[2]^2 + 2.0 * z[1] * z[2] + (z[1]+1.0)^4 + (z[2]-2.0)^4 + z[3]^2 + 2.0 * z[1] * z[3] )
    # flow equations
    @NLconstraint(m, x[1] + x[2] - 1.0 == 0.0)
    @NLconstraint(m, x[1] - x[5] + x[4] == 0.0) #?
    @NLconstraint(m, x[2] - x[6] - x[4] == 0.0)

    # pressure equations
    @NLconstraint(m, x[1]^2 - 1.0 +  y[1] == 0.0)
    @NLconstraint(m, x[2]^2 - 1.0 + y[2] == 0.0)
    @NLconstraint(m, x[4]^2 == y[2] - y[1] ) #??
    @NLconstraint(m, x[5]^2 + y[1] - y[3] == 0.0 )
    @NLconstraint(m, x[6]^2 - y[3] + y[2] == 0.0)

    # stalling constraint
    @constraint(m, z[1]^2 + z[2]^2 + z[3]^2 <= 1.0)

    return m
end

function build_drink3()
    m = Model()
    @variable(m, x[1:3] >= 0.0, start=1.0 )
    @variable(m, y[1:3], start=1.0)
    @variable(m, z[1:3] >= 0.0, start=1.0)

    @NLobjective(m, Min,
        z[1]^2 + 2.0 * z[2]^2 + 2.0 * z[1] * z[2] + (z[1]+0.1)^4 + (z[2]-0.1)^4 + z[3]^2 + 2.0 * z[1] * z[3] )
    # flow equations
    @NLconstraint(m, x[1] + x[2] - 1.0 == 0.0) # supply
    @NLconstraint(m, x[1] + x[3] == 0.5) # demand
    @NLconstraint(m, x[2] - x[3] == 0.5) # demand

    theta = 1.8

    # pressure equations
    @NLconstraint(m, x[1]^theta - 1.0 +  y[1] == 0.0)
    @NLconstraint(m, x[2]^theta - 1.0 + y[2] == 0.0)
    @NLconstraint(m, x[3]^theta == y[1] - y[2] ) #??

    # stalling constraint
    @constraint(m, z[1]^2 + z[2]^2 + z[3]^2 <= 1.0)

    return m
end

function build_drink4()
    m = Model()
    @variable(m, x[1:3] >= 0.0, start=1.0 )
    @variable(m, y[1:3], start=1.0)
    @variable(m, z >= 0.0, start=100.0)

    @NLobjective(m, Min, z^2 )
    # flow equations
    @NLconstraint(m, x[1] + x[2] - 1.0 == 0.0) # supply
    @NLconstraint(m, x[1] + x[3] == 0.5) # demand
    @NLconstraint(m, x[2] - x[3] == 0.5) # demand

    theta = 1.8

    # pressure equations
    @NLconstraint(m, x[1]^theta - 1.0 +  y[1] == 0.0)
    @NLconstraint(m, x[2]^theta - 1.0 + y[2] == 0.0)
    @NLconstraint(m, x[3]^theta == y[1] - y[2] ) #??

    return m
end

function build_drink5()
    m = Model()
    @variable(m, x[1:3] >= 0.0, start=1.0 )
    @variable(m, y[1:3], start=1.0)
    #@variable(m, z >= 0.0, start=100.0)

    @NLobjective(m, Min, 0.0 )
    # flow equations
    @NLconstraint(m, x[1] + x[2] - 1.0 == 0.0) # supply
    @NLconstraint(m, x[1] + x[3] == 0.5) # demand
    @NLconstraint(m, x[2] - x[3] == 0.5) # demand

    theta = 1.8

    # pressure equations
    @NLconstraint(m, x[1]^theta - 1.0 +  y[1] == 0.0)
    @NLconstraint(m, x[2]^theta - 1.0 + y[2] == 0.0)
    @NLconstraint(m, x[3]^theta == y[1] - y[2] ) #??

    return m
end
