# This program is an optimization tool to minimize the cost of
# implementing a submarine grid of dc power lines
# between the Magdalen islands and the mainland

# Define packages used in the program
using JuMP, Clp, Ipopt, AmplNLWriter, CSV, DataFrames

# Define solver used in the optimization problem
# Currently using Couenne as it can handle nonlinear, mixed-integer problems
mod = Model(solver = AmplNLSolver("C:/Users/kevin/Desktop/Design_Project/julia-optimization-tool/Optimization Tool/scipampl_exe/scipampl-6.0.0.win.x86_64.intel.opt.spx2.exe",
["C:/Users/kevin/Desktop/Design_Project//julia-optimization-tool/Optimization Tool/scipampl_exe/scip.set"]))

# Define sets
N = 2 # total number of nodes
T = 100 # largest time value (hour)
L = [0 1; 1 0] # array of possible links (L[n,m] = 1 if there can be links b/w n & m)

# Define parameters
d_nm = [0 241; 241 0] # distances between nodes n & m (km)
a = 1.21e6 # cost of links ($/MW*km)
b = 275e6 # cost of substation ($/MW)
r = 0.009 # resistance on link (ohm/km)
f = 500 # voltage level in (kV), also voltage base
p = 2407 # power capacity of a link (MW)
dem_nt = zeros(Float64, N, T) # power demand @ node n & time t (MW)
for n = 1:N
    for t = 1:T
        if n == 1
            dem_nt[n,t] = 42
        end
        if n == 2
            dem_nt[n,t] = 0
        end
    end
end
lambda_nt = zeros(Float64, N, T) # value of energy @ node & time t ($/MWh)
for n = 1:N
    for t = 1:T
        if n == 1
            lambda_nt[n,t] = 12
        end
        if n == 2
            lambda_nt[n,t] = 7
        end
    end
end
fmax = 1.03 # maximum voltage value for nodes (pu)
fmin = 0.97 # minimum voltage value for nodes (pu)
dr = 0.1 # discount rate
c_n = 500e6*ones(Float64, N, 1) # cost of adding genration @ node n ($)

# Define problem variables
@variables(mod, begin
    x_nm[1:N, 1:N] >= 0, Int # number of parallel lines b/w nodes n & m
    g_nt[1:N, 1:T] >= 0.0 # generation injection @ node n & time t from new generation construction (MW)
    del_nt[1:N, 1:T] >= 0.0 # power injection @ node n & time t from existing generation capacity (MW)
    p_nmt[1:N, 1:N, 1:T] # power flow b/w nodes n & m @ time t (MW)
    u_nt[1:N, 1:T] # voltage @ node n & time t (kV)
    z_n[1:N], Bin # boolean for new generation construction decision
end)

# Define objective function to minimize
@objective(mod, Min, sum(x_nm[n,m]*a*p*d_nm[n,m] for n in 1:N for m in 1:N if n < m) # cost of links
+ sum(x_nm[n,m]*b*p for n in 1:N for m in 1:N if n != m) # cost of substations
+ sum(c_n[n]*z_n[n] for n in 1:N) # cost of generation construction
+ sum(dr^t for t in 1:T)*sum(lambda_nt[n,t]*(g_nt[n,t] - del_nt[n,t]) for n in 1:N for t in 1:T)) # cost of operations

# Add constraints

# Decision to add generation
for n in 1:N
    for t in 1:T
        @constraint(mod, g_nt[n,t] == z_n[n]*g_nt[n,t])
    end
end

# Power balance at a node constraint
for n in 1:N
    for t in 1:T
        @constraint(mod, g_nt[n,t] + del_nt[n,t] == dem_nt[n,t] + sum(p_nmt[n,m,t] for m in 1:N if L[n,m] == 1))
        @constraint(mod, del_nt[1,t] == 0.0) # no ijection from node 1
        @constraint(mod, g_nt[1,t] == 0.0) # no new generation at node 1
    end
end

# Power flow in a link constraint (nonlinear)
for n in 1:N
    for m in 1:N
        for t in 1:T
            if L[n,m] == 1
                @NLconstraint(mod, (p_nmt[n,m,t]) == ((x_nm[n,m]/r)*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t]))
                #@NLconstraint(mod, (p_nmt[n,m,t]) == ((x_nm[n,m]/r_v[v])*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t]))
            end
        end
    end
end

# Power flow bounds constraints
for n in 1:N
    for m in 1:N
        for t in 1:T
            if L[n,m] == 1
                @constraint(mod, (-p*x_nm[n,m]) <= (p_nmt[n,m,t])) # lower bound
                @constraint(mod, (p_nmt[n,m,t]) <= (p*x_nm[n,m])) # upper bound
                @constraint(mod, x_nm[m,n] == x_nm[n,m]) # links b/w n&m = links b/w m&n
            elseif n == m
                @constraint(mod, x_nm[n,m] == 0) # no links between same node
                @constraint(mod, p_nmt[n,m,t] == 0) # no power flow between same node
            end
        end
    end
end

# Voltage bounds constraint
for n in 1:N
    for t in 1:T
        @constraint(mod, (fmin*f) <= (u_nt[n,t]) <= (fmax*f))
    end
end

status = solve(mod) # solve model

# Write outputs to csv files (add file path before file name)
#CSV.write("C:/Users/kevin/Desktop/Design_Project/julia-optimization-tool/Optimization Tool/outputs/links_out.csv", DataFrame(getvalue(x_nm)), append = true)
#CSV.write("gen_sites_out.csv", DataFrame(getvalue(z_n.')), append = true)

# Print results to console
println("x_nm = ", getvalue(x_nm))
#println("g_nt = ", getvalue(g_nt))
#println("del_nt = ", getvalue(del_nt))
#println("p_nmt = ", getvalue(p_nmt))
#println("u_nt = ", getvalue(u_nt))
println("z_n = ", getvalue(z_n.'))
println("Objective value = ", getobjectivevalue(mod))
