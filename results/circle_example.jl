function circle()
    m = Model()
    @variable(m, x >= 0.0, start=1.0)
    @variable(m, y >= 0.0, start=3.0)
    @variable(m, z >= 0.0, start = 2.0)

    @NLobjective(m, Min, x + z^2)
    @NLconstraint(m, x^2 + y^2 <= 1.0)
    @NLconstraint(m, (x-2.0)^2 + y^2 <= 1.0)

    return m
end
