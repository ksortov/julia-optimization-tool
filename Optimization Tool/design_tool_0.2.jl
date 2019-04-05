# This program is an optimization tool to minimize the cost of
# implementing a submarine grid of dc power lines
# between the Magdalen islands and the mainland

# Define packages used in the program
using JuMP, Clp, Ipopt, AmplNLWriter, CSV, DataFrames

# Define solver used in the optimization problem
# Currently using SCIP
mod = Model(solver = AmplNLSolver("C:/Users/kevin/Desktop/Design_Project/julia-optimization-tool/Optimization Tool/scipampl_exe/scipampl-6.0.0.win.x86_64.intel.opt.spx2.exe",
["C:/Users/kevin/Desktop/Design_Project//julia-optimization-tool/Optimization Tool/scipampl_exe/scip.set"]))

inputs = CSV.read("C:/Users/kevin/Desktop/force_wind3.csv") # Read input csv file
init_guess = CSV.read("C:/Users/kevin/Desktop/guesses.csv")
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
lambda_n = inputs[38]
fmax = 1.2 # maximum voltage value for nodes (pu)
fmin = 0.80 # minimum voltage value for nodes (pu)
dr = 0.1 # discount rate
w = 3000 # marginal cost of adding wind genration ($/MW)
c_n = zeros(Float64, N) # variable cost of wind farm construction depending on capacity in MW
for n in 1:N
    c_n[n] = maximum(3000*inputs[n,t+12] for t in 1:T)
end

# Define problem variables
@variables(mod, begin
    x_nm[1:N, 1:N], Bin, Symmetric # number of parallel lines b/w nodes n & m
    0.0 <= g_nt[n in 1:N, t in 1:T] <= inputs[n,t+12] # generation injection @ node n & time t from new generation construction (MW)
    0.0 <= del_nt[n in 1:N, t in 1:T] <= inputs[n,t] # power injection @ node n & time t from existing generation capacity (MW)
    p_nmt[1:N, 1:N, 1:T] # power flow b/w nodes n & m @ time t (MW)
    u_nt[1:N, 1:T] # voltage @ node n & time t (kV)
    z_n[1:N], Bin # boolean for new generation construction decision
    is_sub[1:N], Bin # intermediate variable used to count substations
    sub_num, Int # number of substations to build
    y_v[1:V], Bin # boolean for selected voltage
    alpha_vnm[1:V, 1:N, 1:N], Int # dummy variable for linearization
end)

# Initial guess for the topology x_nm
# for n in 1:N
#     for m in 1:N
#         setValue(x_nm[n,m], init_guess[n,m])
#     end
# end

# Define objective function to minimize
@objective(mod, Min, sum(alpha_vnm[v,n,m]*a_v[v]*d_nm[n,m] for v in 1:V for n in 1:N for m in 1:N if n < m) # cost of links
+ sum(sub_num*y_v[v]*b_v[v] for v in 1:V) # cost of substations
+ sum(c_n[n]*z_n[n] for n in 1:N) # cost of generation construction
+ sum(24*(365/12)*25*lambda_n[n]*(g_nt[n,t] + del_nt[n,t]) for n in 1:N for t in 1:T)) # cost of operations
#sum((1/((1+dr)^t))*(9e6*lambda_n[n]*(g_nt[n,t] + del_nt[n,t])) for n in 1:N for t in 1:T)) # NPV

# Add constraints

# Forced links
#@constraint(mod, x_nm[1,4] == 1)
#@constraint(mod, x_nm[4,1] == 1)
#@constraint(mod, x_nm[1,3] == 1)
#@constraint(mod, x_nm[3,1] == 1)

# Decision to add generation (z_n boolean)
for n in 1:N
    for t in 1:T
        @NLconstraint(mod, g_nt[n,t] == z_n[n]*g_nt[n,t]) # z_n[n] = 1 if g[n,t] > 0 for any t
        @constraint(mod, z_n[n] <= g_nt[n,t]) # ensures z_n[n] = 0 else
            @constraint(mod, g_nt[4,t] >= 0.1)
    end
end

# Number of substations to build (1 per linked node)
for n in 1:N
    for m in 1:N
        @constraint(mod, x_nm[n,m] == is_sub[n]*x_nm[n,m]) # is_sub[n] = 1 if row n in x_nm has at least one nonzero element
    end
end
@constraint(mod, sub_num == sum(is_sub[n] for n in 1:N)) # sum up all elements of is_sub

# Power balance at a node and power injection/ wind generation
for t in 1:T
    for n in 1:N
        @constraint(mod, sum(p_nmt[n,m,t] for m in 1:N if L[n,m] == 1) == g_nt[n,t] + del_nt[n,t] - dem_nt[n,t]) #power balance at a node
    end
end
is

# Power flow in a link constraint (nonlinear)
for n in 1:N
    for m in 1:N
        for t in 1:T
            if n < m
                @NLconstraint(mod, (p_nmt[n,m,t]) == sum((alpha_vnm[v,n,m]/(d_nm[n,m]*r_v[v]))*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t] for v in 1:V)) # power flow equation
                @constraint(mod, p_nmt[n,m,t] == -p_nmt[m,n,t]) # the power from n to m = minus the power for m to n
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
                #@constraint(mod, x_nm[m,n] == x_nm[n,m]) # links b/w n&m = links b/w m&n
            elseif L[n,m] == 0
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
println("sub_num = ", getvalue(sub_num))
println("Objective value = ", getobjectivevalue(mod))

getvalue(x_nm)
getvalue(g_nt)
getvalue(del_nt)
getvalue(p_nmt)
getvalue(u_nt)
getvalue(is_sub)
