# Julia Optimization Tool Utlization guidelines
The following is a list of packages used in the tool (to be added using the Pkg.add("") command):
- JuMP
- Ipopt
- Clp
- AmplNLWriter
The following is an installation guide for Couenne, the utlized solver
- Go to thr opem-source solvers webpage of AMPL Optimization: https://ampl.com/products/solvers/open-source/
- Under "Nonlinear solvers" click the link to download the binaries (not the source code) of the Couenne solver
- In the Julia optimization tool file, in the line: M = Model(solver = AmplNLSolver("path/to/couenne.exe")) replace "path/to/couenne.exe" by the path where you saved couenne.exe