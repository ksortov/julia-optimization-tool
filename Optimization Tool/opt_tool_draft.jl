# This program is an optimization tool to minimize the cost of
# implementing a submarine grid of dc power lines
# between the Magdalen islands and the mainland

# Define packages used in the program
using JuMP, Clp, Ipopt, AmplNLWriter

# Define solver used in the optimization problem
# Currently using Couenne as it can handle nonlinear, mixed-integer problems
#M = Model(solver = IpoptSolver()) # define IpoptSolver as solver (placeholder for now)
#M = Model(solver = AmplNLSolver(Ipopt.amplexe, ["print_level=0 max_cpu_time=30"]))
M = Model(solver = AmplNLSolver("C:/Users/kevin/Desktop/Design Project/couenne-win64/couenne.exe"))

# Define sets
N = 2 # total number of nodes
T = 7 # largest time value (hour)
#V = 4 # number of possible voltage levels
#L = zeros(Int, N, N) # array of possible links (L(n,m) = 1 if there are links b/w n & m)
L = [0 1; 0 0]
#A = 100 # large number used for a constraint on the dummy variable

# Define parameters
#d_nm = ones(Float64, N, N) # distances between nodes n & m (km)
d_nm = [0 241; 0 0]
#a_v = ones(Float64, V, 1) # cost of links ($/MW*km)
a = 1.21e6
#b_v = ones(Float64, V, 1) # cost of substation @ voltage v ($/MW)
b = 275e6
#r_v = ones(Float64, V, 1) # resistance on link for voltage v (ohm/km)
r = 0.009
#f_v = ones(Float64, V, 1) # voltage of v (kV), also voltage base
f = 500
#p_v = ones(Float64, V, 1) # power capacity of a link @ voltage v (MW)
p = 2407
del_nt = zeros(Float64, N, T) # net power injection @ node n & time t (MW)
for n = 1:N
    for t = 1:T
        if n == 1
            del_nt[n,t] = -60
        end
        if n == 2
            del_nt[n,t] = 0
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
fmax = 1.02 # maximum voltage value for nodes (pu)
fmin = 0.98 # minimum voltage value for nodes (pu)
gamma = 0.0175 # discount rate
c_n = zeros(Float64, N, 1) # cost of adding genration @ node n ($)


# Define problem variables
@variable(M, 0 <= x_nm[1:N, 1:N] <= 3, Int) # number of parallel lines b/w nodes n & m
#@variable(M, y_v[1:V], Bin) # boolean for selected voltage
@variable(M, g_nt[1:N, 1:T] >= 0.0) # generation injection @ node n & time t (MW)
@variable(M, p_nmt[1:N, 1:N, 1:T]) # power flow b/w nodes n & m @ time t (MW)
@variable(M, u_nt[1:N, 1:T]) # voltage @ node n & time t (kV)
@variable(M, 0 <= z_n[1:N] <= 1, Int) # boolean for building generation site @ node n
#@variable(M, alpha_vnm[1:V, 1:N, 1:N], Int) # dummy variable for linearization

# Define objective function to minimize
@objective(M, Min, sum(a*p*d_nm[n,m]*x_nm[n,m] for n = 1:N for m = 1:N if n < m) # cost of links
+ sum(b*p*x_nm[n,m] for n = 1:N for m = 1:N if n != m) # cost of substations
+ sum(c_n[n]*z_n[n] for n = 1:N) # cost of generation construction
+ sum(gamma^t for t = 1:T)*sum(lambda_nt[n,t]*(g_nt[n,t] - del_nt[n,t]) for n = 1:N for t = 1:T)) # cost of operations

# Add constraints
# Power balance constraint
for n in 1:N
    for t in 1:T
        @constraint(M, (g_nt[n,t] + del_nt[n,t] - sum(p_nmt[n,m,t] for m = 1:N if (L[n,m] == 1 || L[m,n] == 1))) == 0)
        @constraint(M, g_nt[1,t] == 0)
    end
end


# Power flow in a link constrant
for n in 1:N
    for m in 1:N
        for t in 1:T
            if (L[n,m] == 1 || L[m,n] == 1) #n < m
                #@NLconstraint(M, (p_nmt[n,m,t]) == (sum((alpha_vnm[v,n,m]/r_v[v])*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t] for v = 1:V)))
                @NLconstraint(M, (p_nmt[n,m,t]) == ((x_nm[n,m]/r)*(u_nt[n,t] - u_nt[m,t])*u_nt[n,t]))
                @NLconstraint(M, (p_nmt[m,n,t]) + (p_nmt[n,m,t]) == 0)
            end
        end
    end
end

# Power flow bounds constraints
for n in 1:N
    for m in 1:N
        #for v in 1:V
            for t in 1:T
                #if n < m
                if L[n,m] == 1 #n < m
                    #@constraint(M, (-p_v[v]*x_nm[n,m]) <= (p_nmt[n,m,t]))
                    #@constraint(M, (p_nmt[n,m,t]) <= (p_v[v]*x_nm[n,m]))
                    @constraint(M, (-p*x_nm[n,m]) <= (p_nmt[n,m,t]))
                    @constraint(M, (p_nmt[n,m,t]) <= (p*x_nm[n,m]))
                elseif n == m
                    @constraint(M, x_nm[n,m] == 0)
                end
            end
        #end
    end
end

# Voltage bounds constraint
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

# Print results to console
println("x_nm = ", getvalue(x_nm))
println("g_nt = ", getvalue(g_nt))
println("p_nmt = ", getvalue(p_nmt))
println("u_nt = ", getvalue(u_nt)/f)
println("z_n = ", getvalue(z_n))
println("Objective value = ", getobjectivevalue(M))
