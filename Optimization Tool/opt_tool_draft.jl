using JuMP, Ipopt, Clp

m = Model(solver = IpoptSolver()) # define IpoptSolver as solver (placeholder for now)

# define sets
N = 1 # total number of nodes
T = 1 # largest time value (hour)
V = 1 # number of possible voltage levels
L = Array{Int64}(4, 4) # array of possible links (L(n,m) = 1 if there is a link b/w n % m)
A = 100 # large number used for a constraint on the dummy variable

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
@variable(m, g_nt[1:N, 1:T], Float64) # generation injection @ node n & time t (MW)
@variable(m, p_nmt[1:N, 1:N, 1:T], Float64) # power flow b/w nodes n & m @ time t (MW)
@variable(m, u_nt[1:N, 1:T], Float64) # voltage @ node n & time t (kV)
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
@constraint(m, (g_nt[n,t] - del_nt[n,t] - sum(p_nmt[n,m,t] for m = 1:N) for n = 1:N for t = 1:T) == 0.0)
@constraint(m, (p_nmt[n,m,t] for n = 1:N for m = 1:N for t = 1:T if n !=m) ==
(sum((alpha_vnm[v,n,m]/r_v[v])*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t] for v = 1:V)))
@constraint(m, (-x_nm[n,m]*p_v[v]) <= (p_nmt[n,m,t] for n = 1:N for m = 1:N for v = 1:V for t = 1:t) <= (x_nm[n,m]*p_v[v]))
@constraint(m,(y_v[v]*fmin*f_v[v] for v = 1:V) <= (y_v[v]*fmax*f_v[v]))
@constraint(m, sum(y_v[v] for v = 1:V) == 1)
@constraint(m, (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) == (y_v[v]*x_nm[n,m]))
@constraint(m, 0 <= (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) <= (A*y_v[v]))
@constraint(m, (A*-(1 - y_v[v]) + x_nm[n,m] for n = 1:N for m = 1:N for v = 1:V) <= (alpha_vnm[v,n,m]))
@constraint(m, (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) <= (x_nm[n,m] + A*(1 - y_v[v])))

solve(m)

println("x_nm = ", getvalue(x_nm), " y_v = ", getvalue(y_v),
" g_nt = ", getvalue(g_nt), " p_nmt = ", getvalue(p_nmt),
" u_nt = ", getvalue(u_nt), " z_n = ", getvalue(z_n))
