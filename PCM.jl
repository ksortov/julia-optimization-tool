# Define the packages
using JuMP # used for mathematical programming
using Interact # used for enabling the slider
using Gadfly # used for plotting
using Clp

# Define some input data about the test system
# Maximum power output of generators
const g_max = [0,1000,1000,0,0,0,0,0];
# Minimum power output of generators
const g_min = [0,0,300,0,0,0,0,0];
# Incremental cost of generators
const c_g = [0,50,100,0,0,0,0,0];
# Fixed cost of generators
const c_g0 = [0,1000,10,0,0,0,0,0];
# Incremental cost of wind generators
c_w = [0,0,0,0,0,50,0,0];
# Wind forecast
w_f = [0,0,0,0,0,200,0,0];
#Largest time value
const T=10;

#Resistance
R=0.0366; #ohm/km

#Link lenght
dist=[  0, 223, 329, 235, 196, 117, 257, 198
      223,   0, NaN, NaN, NaN, 275, NaN, NaN
      329, NaN,   0, NaN, NaN, 242, 124, NaN
      235, NaN, NaN,   0, NaN, 216, NaN, 181
      196, NaN, NaN, NaN,   0, 287, NaN, NaN
      117, 275, 242, 216, 287,   0, 167, 115
      257, NaN, 124, NaN, NaN, 167,   0, 183
      198, NaN, NaN, 181, NaN, 115, 183,   0];

#Y matrix
Y = R*dist;

#Number of Nodes
N=8;#not fully sure about this

#Array of possible links
L=[0,1,1,0,0,1,0,0
   1,0,0,0,0,0,0,0
   1,0,0,0,0,0,0,0
   0,0,0,0,0,0,0,0
   0,0,0,0,0,0,0,0
   1,0,0,0,0,0,0,0
   0,0,0,0,0,0,0,0
   0,0,0,0,0,0,0,0];

#Generator Array
#Which nodes can supply power through already existing generation
GB=[0,1,1,0,0,0,0,0];#This is/could essentially be a boolean

#Wind array
#Which nodes can supply Wind
WB=[0,0,0,0,0,1,1,0];#Also essentially a boolean

#voltage levels
v_min=0.95;
v_max=1.05;

#Number of wind turbines
M=1;#this might not work like this for much longer,

#Total demand at every hour
dem=[300, 500, 1500, 900, 1100, 1300, 1500, 1400, 1300, 1200
     0,   0,   0,    0,   0,    0,    0,    0,    0,    0
     0,   0,   0,    0,   0,    0,    0,    0,    0,    0
     0,   0,   0,    0,   0,    0,    0,    0,    0,    0
     0,   0,   0,    0,   0,    0,    0,    0,    0,    0
     0,   0,   0,    0,   0,    0,    0,    0,    0,    0
     0,   0,   0,    0,   0,    0,    0,    0,    0,    0
     0,   0,   0,    0,   0,    0,    0,    0,    0,    0];

g_opt=zeros(N, T);
w_opt=zeros(N, T);
ws_opt=zeros(N, T);
obj=zeros(1,T);

for t in 1:T
    #Define the economic dispatch (ED) model
    ed=Model(solver = ClpSolver())

    # Define decision variables
    @variable(ed, 0 <= g[i=1:N] <= g_max[i]) # power output of generators
    @variable(ed, 0 <= w[i=1:N]  <= w_f[i]) # wind power injection
    @variable(ed, v_nt[i=1:N] <= v_max) # voltage levels


    # Define the objective function
    @objective(ed,Min,sum(c_g[i] * g[i] for i=1:N if GB[i]==1) + sum(c_w[j] * w[j] for j=1:N if WB[j]==1))

    # Define the constraint on the maximum and minimum power output of each generator
    for i in 1:N
        @constraint(ed,  g[i] <= g_max[i]) #maximum
        @constraint(ed,  g[i] >= g_min[i]) #minimum
        # Define the constraint on the wind power injection
        @constraint(ed, w[i] <= w_f[i])
        @constraint(ed, w[i] >= 0)
        #constraints on voltage levels
        @constraint(ed, v_nt[i] <= v_max)
        @constraint(ed, v_nt[i] >= v_min)
    end

    # Define the power balance constraint
    #@constraint(ed, sum(g[i] for i=1:N if GB[i]==1) + sum(w[i] for i=1:N if WB[i]==1) == dem[t])

    #Define power flow/balance
    for i in 1:N
        for k in 1:N
            if L[i,k]==1
                @NLconstraint(ed, sum(v_nt[i] * v_nt[k] * Y[i,k]) == g[i]-dem[i,t])
            end
        end
    end

    # Solve statement
    solve(ed)

    # Solve the economic dispatch problem
    for i=1:N
        if GB[i]==1
            g_opt[i,t]=getvalue(g[i])
        elseif WB[i]==1
            w_opt[i,t]=getvalue(w[i])
            ws_opt[i,t]=w_f[i]-getvalue(w[i])
        end
    end
    obj[t]=getobjectivevalue(ed)
end

T_cost=sum(obj[t] for t=1:T)

println("Dispatch of Generators: ", g_opt, " MW");
println("Dispatch of Wind: ", w_opt, " MW");
println("Wind spillage: ", ws_opt, " MW");
println("\n");
println("Hourly cost: ", obj, " \$");
println("Total cost: ", T_cost, "\$");
g_opt
w_opt
obj
T_cost