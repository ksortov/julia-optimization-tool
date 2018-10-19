using JuMP
using Ipopt

m = Model(solver = IpoptSolver()) # define IpoptSolver as solver (placeholder for now)

# define sets
N = 4 # total number of nodes
T = 1 # largest time value (hour)
V = 1 # number of possible voltage levels
L = Array{Int64}(4, 4) # array of possible links (L(n,m) = 1 if there is a link b/w n % m)

# define parameters
d_nm = Array{Float64}(N, N) # distances between nodes n & m (km)
a_v = Array{Float64}(1, V) # cost of links ($/MW*km)
b_v = Array{Float64}(1, V) # cost of substation @ voltage v ($/MW)
r_v = Array{Float64}(1, V) # resistance on link for voltage v (ohm/km)
f_v = Array{Float64}(1, V) # voltage of v (kV)
p_v = Array{Float64}(1, V) # power capacity of a link @ voltage v (MW)
del_nt = Array{Float64}(N, T) # net power injection @ node n & time t (MW)
lambda_nt = Array{Float64}(N, T) # value of energy @ node & time t ($/MWh)
fmax = 1 # maximum voltage value for nodes (pu)
fmin = 1 # minimum voltage value for nodes (pu)
gamma = 1 # discount rate
c_n = Array{Float64}(1, N) # cost of adding genration @ node n ($)

# define problem variables
@variable(m, x_nm[1:N, 1:N], Int) # number of parallel lines b/w nodes n & m
@variable(m, y_v[1:V], Bin) # boolean for selected voltage
@variable(m, g_nt[1:N, 1:T]) # generation injection @ node n & time t (MW)
@variable(m, p_nmt[1:N, 1:N, 1:T]) # power flow b/w nodes n & m @ time t (MW)
@variable(m, u_nt[1:N, 1:T]) # voltage @ node n & time t (kV)
@variable(m, z_n[1:N], Bin) # boolean for building generation site @ node n
@variable(m, alpha_vnm[1:V, 1:N, 1:N], Int) # dummy variable for linearization

# define equations to be used in the objective function
# capital equations = sum of capital costs for links, substations, generation construction

# cost for links
@expression(m, links, sum(a_v[v]*p_v[v]*d_nm[n,m]*y_v[v]*x_nm[n,m] for n = 1:N for m = 1:N for v = 1:V if n < m))

# cost for substations
@expression(m, subs, sum(b_v[v]*p_v[v]*y_v[v] for v = 1:V)*sum(x_nm[n,m] for n = 1:N for m = 1:N if n != m))

# cost for generation construction
@expression(m, gens, sum(c_n[n]*z_n[n] for n = 1:N))

# operations cost equation
@expression(m, ops, sum(gamma^t for t = 1:T)*sum(lambda_nt[n,t]*(g_nt[n,t] - del_nt[n,t]) for n = 1:N for t = 1:T))

#@NLobjective(m, Min, links + subs + gens + ops) # define function to optimize

# adding constraints
@constraint(m, (g_nt - del_nt - sum(p_nmt[n,m,t] for m = 1:N) for n = 1:N for t = 1:T) == 0.0)
@constraint(m,)
@constraint(m,)
@constraint(m,)
@constraint(m, sum(y_v[v] for v = 1:V) == 1)
@constraint(m,)

solve(m)

println("x_nm = ", getvalue(x_nm), " y_v = ", getvalue(y_v),
" g_nt = ", getvalue(g_nt), " p_nmt = ", getvalue(p_nmt),
" u_nt = ", getvalue(u_nt), " z_n = ", getvalue(z_n))
