function simple_comp1()
    m = Model()
    @variable(m, x >= 0.0)
    @variable(m, y >= 0.0)
    @variable(m, l[1:3] >= 0.0)
    @variable(m, z[1:3] >= 0.0)

    @NLobjective(m, Min, (x - 5)^2 + (2*y + 1)^2)

    @NLconstraint(m, 2*(y-1) - 1.5*x + l[1] - l[2]*0.5 + l[3] == 0)

    @constraint(m, 3*x - y - 3 == z[1])
    @constraint(m, - x + 0.5*y + 4 == z[2])
    @constraint(m, - x - y + 7 == z[3])
    for i = 1:3
        @NLconstraint(m, z[i] * l[i] <= 0.0)
    end

    return m
end

function simple_comp2()
    m = Model()
    @variable(m, x >= 0.0)
    @variable(m, 0.5 >= y >= 0.0)

    @NLobjective(m, Min, 3 * x + y)

    @constraint(m, x + y == 1.0)
    @NLconstraint(m, x * y <= 0.0)

    return m
end
