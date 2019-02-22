# This program is an optimization tool to minimize the cost of
# implementing a submarine grid of dc power lines
# between the Magdalen islands and the mainland

# Define packages used in the program
using JuMP, Clp, Ipopt, AmplNLWriter, CSV, DataFrames

# Define solver used in the optimization problem
# Currently using Couenne as it can handle nonlinear, mixed-integer problems
mod = Model(solver = AmplNLSolver("C:/Users/kevin/Desktop/Design_Project/julia-optimization-tool/Optimization Tool/scipampl_exe/scipampl-6.0.0.win.x86_64.intel.opt.spx2.exe",
["C:/Users/kevin/Desktop/Design_Project//julia-optimization-tool/Optimization Tool/scipampl_exe/scip.set"]))

inputs = CSV.read("C:/Users/kevin/Desktop/scenarios/scenario6.csv") # Read input csv file
node_num = length(inputs[1]) # number of nodes considered in given scenario

# Define sets
V = 2 # total number of potential voltage levels
N = node_num # total number of nodes
T = 12 # largest time value (hour)
A = 1000 # large number used for dummy variable constraints
L = inputs[47:47+node_num-1] # array of possible links (L[n,m] = 1 if there can be links b/w n & m)

# Define parameters
d_nm = inputs[47+node_num:47+2*node_num-1] # distances between nodes n & m (km)
a_v = inputs[40] # cost of links ($/km)
b_v = inputs[41] # cost of substation ($/station)
r_v = inputs[42] # resistance on link (ohm/km)
f_v = inputs[39] # voltage level in (kV), also voltage base
p_v = inputs[43] # power capacity of a link (MW)
dem_nt = inputs[25:36] # power demand @ node n & time t (MW)
lambda_nt = zeros(Float64, N, T) # value of energy @ node & time t ($/MWh)
for n = 1:N
    for t = 1:T
            lambda_nt[n,t] = inputs[n,38]
    end
end
fmax = 1.05 # maximum voltage value for nodes (pu)
fmin = 0.95 # minimum voltage value for nodes (pu)
dr = 0.1 # discount rate
#c_n = [0 0 30e3] # cost of adding wind generation at node n ($)
w = 3000 # marginal cost of adding wind genration ($/MW)

# Define problem variables
@variables(mod, begin
    x_nm[1:N, 1:N] >= 0, Int # number of parallel lines b/w nodes n & m
    g_nt[1:N, 1:T] >= 0.0 # generation injection @ node n & time t from new generation construction (MW)
    del_nt[1:N, 1:T] >= 0.0 # power injection @ node n & time t from existing generation capacity (MW)
    p_nmt[1:N, 1:N, 1:T] # power flow b/w nodes n & m @ time t (MW)
    u_nt[1:N, 1:T] # voltage @ node n & time t (kV)
    c_n[1:N] # total cost of adding wind generation at node n ($)
    z_n[1:N], Bin # boolean for new generation construction decision
    y_v[1:V], Bin # boolean for selected voltage
    alpha_vnm[1:V, 1:N, 1:N], Int # dummy variable for linearization
end)

# Define objective function to minimize
@objective(mod, Min, sum(alpha_vnm[v,n,m]*a_v[v]*d_nm[n,m] for v in 1:V for n in 1:N for m in 1:N if L[n,m] == 1) # cost of links
+ sum(alpha_vnm[v,n,m]*b_v[v] for v in 1:V for n in 1:N for m in 1:N if n != m) # cost of substations
+ sum(c_n[n]*z_n[n] for n in 1:N) # cost of generation construction
+ sum(dr^t for t in 1:T)*sum(lambda_nt[n,t]*(g_nt[n,t] + del_nt[n,t]) for n in 1:N for t in 1:T)) # cost of operations

# Add constraints

for n in 1:N
    @constraint(mod, c_n[n] == maximum(w*g_nt[n,]))
end

# Decision to add generation (z_n boolean)
for n in 1:N
    for t in 1:T
        @constraint(mod, g_nt[n,t] == z_n[n]*g_nt[n,t])
        @constraint(mod, z_n[n] <= g_nt[n,t])
    end
end

# Power balance at a node and power injection/ wind generation
for t in 1:T
    for n in 1:N
        @constraint(mod, g_nt[n,t] <= inputs[n,t+12]) # maximum amount of wind generation @ node n and time t
        @constraint(mod, g_nt[n,t] + del_nt[n,t] == dem_nt[n,t] + sum(p_nmt[n,m,t] for m in 1:N if n != m)) #power balance at a node
        if n != 2
            @constraint(mod, del_nt[n,t] == inputs[n,t]) # injection from node n (node 2 is slack)
        end
    end
    @constraint(mod, sum(p_nmt[n,m,t] for n in 1:N for m in 1:N) == 0) # power balance of the system
end

# Power flow in a link constraint (nonlinear)
for n in 1:N
    for m in 1:N
        for t in 1:T
            if L[n,m] == 1
                @NLconstraint(mod, (p_nmt[n,m,t]) == sum((alpha_vnm[v,n,m]/r_v[v])*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t] for v in 1:V))
            end
        end
    end
end

# Power flow bounds constraints
for n in 1:N
    for m in 1:N
        for t in 1:T
            if L[n,m] == 1
                @constraint(mod, (p_nmt[n,m,t]) >= sum(-p_v[v]*x_nm[n,m] for v in 1:V)) # lower bound
                @constraint(mod, (p_nmt[n,m,t]) <= sum(p_v[v]*x_nm[n,m] for v in 1:V)) # upper bound
                @constraint(mod, x_nm[m,n] == x_nm[n,m]) # links b/w n&m = links b/w m&n
            elseif n == m || L[n,m] == 0
                @constraint(mod, x_nm[n,m] == 0) # no links between same node
                @constraint(mod, p_nmt[n,m,t] == 0) # no power flow between same node
            end
        end
    end
end

# Voltage bounds constraint
for n in 1:N
    for t in 1:T
        @constraint(mod, (u_nt[n,t]) >= sum(fmin*f_v[v]*y_v[v] for v in 1:V)) # lower bound
        @constraint(mod, (u_nt[n,t]) <= sum(fmax*f_v[v]*y_v[v] for v in 1:V)) # upper bound
    end
end

# Logical constraint on y_v (only one voltage level is selected)
@constraint(mod, (sum(y_v[v] for v in 1:V)) == 1)

# Constarints on linearization variable alpha
# which is equal to x_nm when the boolean y_v is 1
for v in 1:V
    for n in 1:N
        for m in 1:N
            @NLconstraint(mod, (alpha_vnm[v,n,m]) == (y_v[v]*x_nm[n,m])) # alpha is the product of links and selected voltage boolean
            @constraint(mod, (alpha_vnm[v,n,m]) >= 0) # minimum value is 0
            @constraint(mod, (alpha_vnm[v,n,m]) <= (A*y_v[v])) # maximum value is a large number
            @constraint(mod, (alpha_vnm[v,n,m]) >= (A*-(1 - y_v[v]) + x_nm[n,m])) # these two last constranits ensure that alpha
            @constraint(mod, (alpha_vnm[v,n,m]) <= (x_nm[n,m] + A*(1 - y_v[v])))  # actually takes on the value of x_nm (links)
        end
    end
end

# N-1 criterion
#for n in 1:N
#    @constraint(mod, sum(x_nm[1,m] for m in 2:N) >= 2)
#end

status = solve(mod) # solve model

# Write outputs to csv files (add file path before file name)
#CSV.write("C:/Users/kevin/Desktop/Design_Project/julia-optimization-tool/Optimization Tool/outputs/links_out.csv", DataFrame(getvalue(x_nm)), append = true)
#CSV.write("C:/Users/kevin/Desktop/links_out.csv", DataFrame(getvalue(x_nm)), append = true)

# Print results to console
println("y_v = ", getvalue(y_v))
println("x_nm = ", getvalue(x_nm))
println("g_nt = ", getvalue(g_nt))
println("del_nt = ", getvalue(del_nt))
println("p_nmt = ", getvalue(p_nmt))
println("u_nt = ", getvalue(u_nt))
println("z_n = ", getvalue(z_n.'))
println("Objective value = ", getobjectivevalue(mod))
