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
# Total demand
const d = 300;
# Wind forecast
const w_f = 200;


#Define the economic dispatch (ED) model
ed=Model(solver = ClpSolver())

# Define decision variables
@variable(ed, 0 <= g[i=1:2] <= g_max[i]) # power output of generators
@variable(ed, 0 <= w  <= w_f ) # wind power injection

# Define the objective function
@objective(ed,Min,sum(c_g[i] * g[i] for i=1:2)+ c_w * w)

# Define the constraint on the maximum and minimum power output of each generator
for i in 1:2
    @constraint(ed,  g[i] <= g_max[i]) #maximum
    @constraint(ed,  g[i] >= g_min[i]) #minimum
end

# Define the constraint on the wind power injection
@constraint(ed, w <= w_f)

# Define the power balance constraint
@constraint(ed, sum(g[i] for i=1:2) + w == d)

# Solve statement
solve(ed)



# Solve the economic dispatch problem
g_opt=getValue(g)
w_opt=getValue(w)
ws_opt=w_f-getValue(w)
obj=getObjectiveValue(ed)

#(g_opt,w_opt,ws_opt,obj) = solve_ed(g_max, g_min, c_g, c_w, d, w_f);

println("\n")
println("Dispatch of Generators: ", g_opt[i=1:2], " MW")
println("Dispatch of Wind: ", w_opt, " MW")
println("Wind spillage: ", w_f-w_opt, " MW")
println("\n")
println("Total cost: ", obj, "\$")


using DataFrames
using CSV

T=24;
N=8;

g_opt=ones(N,T);
w_opt=zeros(N,T);
v_opt=50*ones(N,T);
obj=2*ones(1,T);

out_file=[  "Hydro generation" fill(-, 1, T-1)
            g_opt
            "Wind generation" fill(-, 1, T-1)
            w_opt
            "Voltage levels" fill(-, 1, T-1)
            v_opt
            "Hourly cost" fill(-, 1, T-1)
            obj]

# Write outputs to csv files (add file path before file name)
#CSV.write("C:/Users/kevin/Desktop/Design_Project/julia-optimization-tool/Optimization Tool/outputs/links_out.csv", DataFrame(getvalue(x_nm)), append = true)
CSV.write("/Users/Antoine/Documents/pcm_out.csv", DataFrame(out_file), append = true)
