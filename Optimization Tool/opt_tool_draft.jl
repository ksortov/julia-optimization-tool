using JuMP
using Ipopt

m = Model(solver = IpoptSolver()) # define IpoptSolver as solver

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
c_n = 1 # cost of adding genration @ node n ($)

# define problem variables
@variable(m, x_nm[1:N, 1:N], Int) # number of parallel lines b/w nodes n & m
@variable(m, y_v[1:V], Bin) # boolean for selected voltage
@variable(m, g_nt[1:N, 1:T], start = 0.0) # generation injection @ node n & time t (MW)
@variable(m, p_nmt[1:N, 1:N, 1:T], start = 0.0) # power flow b/w nodes n & m @ time t (MW)
@variable(m, u_nt[1:N, 1:T], start = 0.0) # voltage @ node n & time t (kV)
@variable(m, z_n[1:N], Bin) # boolean for building generation site @ node n
@variable(m, alpha_vnm, Int) # dummy variable for linearization

@NLobjective(m, Min, ) # define function to optimize

# adding constraints
@constraint(m,)
@constraint(m,)
@constraint(m,)
@constraint(m,)
@constraint(m, sum{y_v[i], i = V} == 1)
@constraint(m,)

solve(m)

println("x_nm = ", getvalue(x_nm), " y_v = ", getvalue(y_v),
" g_nt = ", getvalue(g_nt), " p_nmt = ", getvalue(p_nmt),
" u_nt = ", getvalue(u_nt), " z_n = ", getvalue(z_n))
