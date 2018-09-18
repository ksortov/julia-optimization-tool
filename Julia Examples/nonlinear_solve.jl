# define packages to be used
using JuMP
using Ipopt

m = Model(solver = IpoptSolver()) # define IpoptSolver as solver

# define problem variables
@variable(m, x, start = 0.0)
@variable(m, y, start = 0.0)

@NLobjective(m, Min, (1-x)^2 + 100(y-x^2)^2) # define function to optimize

# solve and print solution
solve(m)
println("x = ", getvalue(x), " y = ", getvalue(y))

# adding a (linear) constraint
@constraint(m, x + y == 10)
solve(m)
println("x = ", getvalue(x), " y = ", getvalue(y))
