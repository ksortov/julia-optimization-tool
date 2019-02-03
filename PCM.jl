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

#Number of wind turbines
M=1;

#Total demand at every hour
dem=[300, 500, 700, 900, 1100, 1300, 1500, 1400, 1300, 1200];

for t in 1:2

    de=dem[t];

    #Define the economic dispatch (ED) model
    ed=Model(solver = ClpSolver())

    # Define decision variables
    @variable(ed, 0 <= g[1:N] <= g_max[i]) # power output of generators
    @variable(ed, 0 <= w[1:M]  <= w_f ) # wind power injection

    # Define the objective function
    @objective(ed,Min,sum(c_g[i] * g[i] for i=1:N) + sum(c_w[j] * w[j] for j=1:M))

    # Define the constraint on the maximum and minimum power output of each generator
    for i in 1:N
        @constraint(ed,  g[i] <= g_max[i]) #maximum
        @constraint(ed,  g[i] >= g_min[i]) #minimum
    end

    for j in 1:M
        # Define the constraint on the wind power injection
        @constraint(ed, w[j] <= w_f)
        @constraint(ed, w[j] >= 0)
    end

    # Define the power balance constraint
    @constraint(ed, sum(g[i] for i=1:N) + sum(w[j] for j=1:M) == de)

    # Solve statement
    solve(ed)

    # Solve the economic dispatch problem
    for i=1:N
        g_opt[i,t]=getvalue(g[i])
    end
    for j=1:M
        w_opt[j,t]=getvalue(w[j])
        #ws_opt[j,t]=w_f-getvalue(w[j])
    end
    #obj[t]=getobjectivevalue(ed)
end

#T_cost=sum(obj[t] for t=1:T)

println("Dispatch of Generators: ", g_opt, " MW")
println("Dispatch of Wind: ", w_opt, " MW")
println("Wind spillage: ", ws_opt, " MW")
println("\n")
println("Hourly cost: ", obj, " \$")
println("Total cost: ", T_cost, "\$")
