# Julia Optimization Tool Utilization Guidelines
The following is a list of packages used in the tool (to be added using the Pkg.add("") command in the Julia console):
- JuMP
- Ipopt
- Clp
- AmplNLWriter
- CSV
- DataFrames

The following is a guide on how to use SCIP, the utilized solver:
- The solver executable file for windows is located in the scipampl_exe directory
- In the Julia optimization tool file, in the line: mod = Model(solver = AmplNLSolver("path/to/scip.exe")) replace "path/to/scip.exe" by your local file path to the couenne executable
