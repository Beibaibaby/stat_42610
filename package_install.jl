using Pkg
#Pkg.add("Distributions")
#Pkg.add("ProgressMeter")
#Pkg.add("Plots")
#Pkg.add("JSON3")
#Pkg.add("JSON3")
#Pkg.add("Measures")
#Pkg.add("JLD2")
#Pkg.add("SharedArrays")
#Pkg.add("Distributed")
#Pkg.update()
#Pkg.build("GR")
#Pkg.add("FilePathsBase")
#Pkg.add("FileIO")
#Pkg.add("DSP")
#Pkg.add("Optim")
#Pkg.add("DifferentialEquations")
#
#Pkg.add("Conda") 
#using Conda
#Pkg.add("SymPy")
#Conda.update()
using Conda
Conda.add("sympy")

