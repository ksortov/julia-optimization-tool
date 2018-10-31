using JuMP, Ipopt, Clp, NLopt

M = Model(solver = IpoptSolver()) # define IpoptSolver as solver (placeholder for now)
#m = Model(solver=NLoptSolver(algorithm=:LD_MMA))

# define sets
N = 2 # total number of nodes
T = 3 # largest time value (hour)
V = 4 # number of possible voltage levels
L = ones(Int, 4, 4) # array of possible links (L(n,m) = 1 if there is a link b/w n & m)
A = 100 # large number used for a constraint on the dummy variable

# define parameters
d_nm = ones(Float64, N, N) # distances between nodes n & m (km)
a_v = ones(Float64, V, 1) # cost of links ($/MW*km)
b_v = ones(Float64, V, 1) # cost of substation @ voltage v ($/MW)
r_v = ones(Float64, V, 1) # resistance on link for voltage v (ohm/km)
f_v = ones(Float64, V, 1) # voltage of v (kV)
p_v = ones(Float64, V, 1) # power capacity of a link @ voltage v (MW)
del_nt = zeros(Float64, N, T) # net power injection @ node n & time t (MW)
lambda_nt = ones(Float64, N, T) # value of energy @ node & time t ($/MWh)
fmax = 1 # maximum voltage value for nodes (pu)
fmin = 1 # minimum voltage value for nodes (pu)
gamma = 1 # discount rate
c_n = ones(Float64, N, 1) # cost of adding genration @ node n ($)


# define problem variables
@variable(M, x_nm[1:N, 1:N], Int) # number of parallel lines b/w nodes n & m
@variable(M, y_v[1:V], Bin) # boolean for selected voltage
@variable(M, g_nt[1:N, 1:T]) # generation injection @ node n & time t (MW)
@variable(M, p_nmt[1:N, 1:N, 1:T]) # power flow b/w nodes n & m @ time t (MW)
@variable(M, u_nt[1:N, 1:T]) # voltage @ node n & time t (kV)
@variable(M, z_n[1:N], Bin) # boolean for building generation site @ node n
@variable(M, alpha_vnm[1:V, 1:N, 1:N], Int) # dummy variable for linearization


# define equations to be used in the objective function
# capital equations = sum of capital costs for links, substations, generation construction

# cost for links
@expression(M, links, sum(a_v[v]*p_v[v]*d_nm[n,m]*alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V if n < m))

# cost for substations
@expression(M, subs, sum(b_v[v]*p_v[v]*y_v[v] for v = 1:V)*sum(x_nm[n,m] for n = 1:N for m = 1:N if n != m))

# cost for generation construction
@expression(M, gens, sum(c_n[n]*z_n[n] for n = 1:N))

# operations cost equation
@expression(M, ops, sum(gamma^t for t = 1:T)*sum(lambda_nt[n,t]*(g_nt[n,t] - del_nt[n,t]) for n = 1:N for t = 1:T))

@objective(M, Min, links + subs + gens + ops) # define objective function

# adding constraints
for n in 1:N
    for t in 1:T
        @constraint(M, (g_nt[n,t] - del_nt[n,t] - sum(p_nmt[n,m,t] for m = 1:N)) == 0)
    end
end

for n in 1:N
    for m in 1:N
        for t in 1:T
            if n !=m
                @NLconstraint(M, (p_nmt[n,m,t]) == (sum((alpha_vnm[v,n,m]/r_v[v])*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t] for v = 1:V)))
            end
        end
    end
end

for n in 1:N
    for m in 1:N
        for v in 1:V
            for t in 1:T
                if n < m && n < v && n < t
                    @constraint(M, (-p_v[v]*x_nm[n,m]) <= (p_nmt[n,m,t]) <= (p_v[v]*x_nm[n,m]))
                end
            end
        end
    end
end
# @constraint(m,(y_v[v]*fmin*f_v[v] for v = 1:V) <= (y_v[v]*fmax*f_v[v]))
# @constraint(m, sum(y_v[v] for v = 1:V) == 1)
# @constraint(m, (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) == (y_v[v]*x_nm[n,m]))
# @constraint(m, 0 <= (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) <= (A*y_v[v]))
# @constraint(m, (A*-(1 - y_v[v]) + x_nm[n,m] for n = 1:N for m = 1:N for v = 1:V) <= (alpha_vnm[v,n,m]))
# @constraint(m, (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) <= (x_nm[n,m] + A*(1 - y_v[v])))

#
# solve(m)
#
# println("alpha_nm = ", getvalue(x_nm), " y_v = ", getvalue(y_v),
# " g_nt = ", getvalue(g_nt), " p_nmt = ", getvalue(p_nmt),
# " u_nt = ", getvalue(u_nt), " z_n = ", getvalue(z_n))
