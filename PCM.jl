# Define the packages
using JuMP # used for mathematical programming
using Interact # used for enabling the slider
using Gadfly # used for plotting
using Clp

# Define some input data about the test system
# Maximum power output of generators
const g_max = [1000,1000];
# Minimum power output of generators
const g_min = [0,300];
# Incremental cost of generators
const c_g = [50,100];
# Fixed cost of generators
const c_g0 = [1000,0]
# Incremental cost of wind generators
const c_w = 50;
# Wind forecast
const w_f = 200;
#Largest time value
const T=10;

#Number of Generators
N=2;

#Total demand at every hour
dem=[300, 500, 700, 900, 1100, 1300, 1500, 1400, 1300, 1200];

#Define the economic dispatch (ED) model
ed=Model(solver = ClpSolver())


# Define decision variables
@variable(ed, 0 <= g[1:N, 1:T] <= g_max[i]) # power output of generators
@variable(ed, 0 <= w[1:T]  <= w_f ) # wind power injection


# Define the objective function
for t in 1:T
    @objective(ed,Min,sum(c_g[i] * g[i,t] for i=1:N)+ c_w * w[t])
end

# Define the constraint on the maximum and minimum power output of each generator
for i in 1:N
    for t in 1:T
        @constraint(ed,  g[i,t] <= g_max[i]) #maximum
        @constraint(ed,  g[i,t] >= g_min[i]) #minimum
    end
end

# Define the constraint on the wind power injection
for t in 1:T
    @constraint(ed, w[t] <= w_f)
end

# Define the power balance constraint
for t=1:T
    @constraint(ed, sum(g[i,t] for i=1:N) + w[t] == dem[t])
end

# Solve statement
solve(ed)



# Solve the economic dispatch problem
g_opt=getValue(g)
w_opt=getValue(w)
ws_opt=w_f-getValue(w)
obj=getObjectiveValue(ed)

#(g_opt,w_opt,ws_opt,obj) = solve_ed(g_max, g_min, c_g, c_w, d, w_f);

println("\n")
println("Dispatch of Generators: ", g_opt, " MW")
println("Dispatch of Wind: ", w_opt, " MW")
println("Wind spillage: ", w_f-w_opt, " MW")
println("\n")
println("Total cost: ", obj, "\$")
