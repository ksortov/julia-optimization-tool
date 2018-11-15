# Julia Optimization Tool Utlization Guidelines
The following is a list of packages used in the tool (to be added using the Pkg.add("") command in the Julia console):
- JuMP
- Ipopt
- Clp
- AmplNLWriter

The following is a guide on how to use Couenne, the utilized solver:
- The solver executable is located in the couenne-win64 directory
- In the Julia optimization tool file, in the line: M = Model(solver = AmplNLSolver("path/to/couenne.exe")) replace "path/to/couenne.exe" by your local file path to couenne.exe
