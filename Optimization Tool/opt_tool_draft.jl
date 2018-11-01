using JuMP, Ipopt, Clp, AmplNLWriter, CoinOptServices

#M = Model(solver = IpoptSolver()) # define IpoptSolver as solver (placeholder for now)
#M = Model(solver = AmplNLSolver(Ipopt.amplexe, ["print_level=0 max_cpu_time=30"]))
M = Model(solver = AmplNLSolver("C:/Users/kevin/Desktop/Design Project/couenne-win64/couenne.exe"))

# define sets
N = 9 # total number of nodes
T = 5 # largest time value (hour)
V = 4 # number of possible voltage levels
L = ones(Int, 4, 4) # array of possible links (L(n,m) = 1 if there is a link b/w n & m)
A = 100 # large number used for a constraint on the dummy variable

# define parameters
#d_nm = ones(Float64, N, N) # distances between nodes n & m (km)
d_nm = [0 112 0 143 175 0 0 0 0
0 0 89 0 140 0 0 0 0
0 0 0 0 161 196 0 0 0
0 0 0 0 120 0 194 0 0
0 0 0 0 0 147 288 177 269
0 0 0 0 0 0 0 0 170
0 0 0 0 0 0 0 241 0
0 0 0 0 0 0 0 0 218
0 0 0 0 0 0 0 0 0]
#a_v = ones(Float64, V, 1) # cost of links ($/MW*km)
a = 1e6
#b_v = ones(Float64, V, 1) # cost of substation @ voltage v ($/MW)
b = 1
#r_v = ones(Float64, V, 1) # resistance on link for voltage v (ohm/km)
r = 0.012
#f_v = ones(Float64, V, 1) # voltage of v (kV)
f = 320
#p_v = ones(Float64, V, 1) # power capacity of a link @ voltage v (MW)
p = 500
del_nt = ones(Float64, N, T) # net power injection @ node n & time t (MW)
lambda_nt = ones(Float64, N, T) # value of energy @ node & time t ($/MWh)
fmax = 1.05 # maximum voltage value for nodes (pu)
fmin = 0.95 # minimum voltage value for nodes (pu)
gamma = 1 # discount rate
c_n = ones(Float64, N, 1) # cost of adding genration @ node n ($)


# define problem variables
@variable(M, x_nm[1:N, 1:N], Int) # number of parallel lines b/w nodes n & m
#@variable(M, y_v[1:V], Bin) # boolean for selected voltage
@variable(M, g_nt[1:N, 1:T]) # generation injection @ node n & time t (MW)
@variable(M, p_nmt[1:N, 1:N, 1:T]) # power flow b/w nodes n & m @ time t (MW)
@variable(M, u_nt[1:N, 1:T]) # voltage @ node n & time t (kV)
@variable(M, 0 <= z_n[1:N] <= 1, Int) # boolean for building generation site @ node n
#@variable(M, alpha_vnm[1:V, 1:N, 1:N], Int) # dummy variable for linearization


# define equations to be used in the objective function
# capital equations = sum of capital costs for links, substations, generation construction

# cost for links
#@expression(M, links, sum(a_v[v]*p_v[v]*d_nm[n,m]*alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V if n < m))
@expression(M, links, sum(a*p*d_nm[n,m]*x_nm[n,m] for n = 1:N for m = 1:N if n < m))

# cost for substations
#@expression(M, subs, sum(b_v[v]*p_v[v]*y_v[v] for v = 1:V)*sum(x_nm[n,m] for n = 1:N for m = 1:N if n != m))
@expression(M, subs, sum(b*p*x_nm[n,m] for n = 1:N for m = 1:N if n != m))

# cost for generation construction
@expression(M, gens, sum(c_n[n]*z_n[n] for n = 1:N))

# operations cost equation
@expression(M, ops, sum(gamma^t for t = 1:T)*sum(lambda_nt[n,t]*(g_nt[n,t] - del_nt[n,t]) for n = 1:N for t = 1:T))

# sum of all cost functions
@expression(M, objective, links + subs + gens + ops)

@objective(M, Min, objective) # define objective function

# adding constraints
# lower bound on x_nm
for n in 1:N
    for m in 1:N
        @constraint(M, x_nm[n,m] >= 0) # lower bound on objective value
    end
end

# power balance cnstraint
for n in 1:N
    for t in 1:T
        @constraint(M, (g_nt[n,t] - del_nt[n,t] - sum(p_nmt[n,m,t] for m = 1:N)) == 0)
    end
end


# power flow in a link constrant
for n in 1:N
    for m in 1:N
        for t in 1:T
            if n !=m
                #@NLconstraint(M, (p_nmt[n,m,t]) == (sum((alpha_vnm[v,n,m]/r_v[v])*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t] for v = 1:V)))
                @NLconstraint(M, (p_nmt[n,m,t]) == ((x_nm[n,m]/r)*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t]))
            end
        end
    end
end

# power flow bounds constraints
for n in 1:N
    for m in 1:N
        #for v in 1:V
            for t in 1:T
                #if n < m && n < v && n < t
                if n < m && n < t
                    #@constraint(M, (-p_v[v]*x_nm[n,m]) <= (p_nmt[n,m,t]))
                    #@constraint(M, (p_nmt[n,m,t]) <= (p_v[v]*x_nm[n,m]))
                    @constraint(M, (-p*x_nm[n,m]) <= (p_nmt[n,m,t]))
                    @constraint(M, (p_nmt[n,m,t]) <= (p*x_nm[n,m]))
                end
            end
        #end
    end
end

# voltage bounds constraint
for n in 1:N
    for t in 1:T
        #for v in 1:V
        #@constraint(M,(y_v[v]*fmin*f_v[v]) <= (y_v[v]*fmax*f_v[v]))
        @constraint(M, (fmin*f) <= (u_nt[n,t]) <= (fmax*f))
        #end
    end
end
# logical constraint on y_v
#@constraint(M, sum(y_v[v] for v = 1:V) == 1)

# constarints on dummy alpha variable (not used for now)
# @constraint(M, (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) == (y_v[v]*x_nm[n,m]))
# @constraint(m, 0 <= (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) <= (A*y_v[v]))
# @constraint(m, (A*-(1 - y_v[v]) + x_nm[n,m] for n = 1:N for m = 1:N for v = 1:V) <= (alpha_vnm[v,n,m]))
# @constraint(m, (alpha_vnm[v,n,m] for n = 1:N for m = 1:N for v = 1:V) <= (x_nm[n,m] + A*(1 - y_v[v])))


solve(M) # solve model

# print results to console
println("x_nm = ", getvalue(x_nm))
println("g_nt = ", getvalue(g_nt))
println("p_nmt = ", getvalue(p_nmt))
println(" u_nt = ", getvalue(u_nt))
println(" z_n = ", getvalue(z_n))
println("Objective value: ", getobjectivevalue(M))
