# To include this file:
# include("Box/FLOW Lab/Modules MCA/NonlinearLiftingLine/src/NLL Test.jl")

#Pkg.develop(PackageSpec(path = "/Users/markanderson/Box/FLOW Lab/Modules MCA/makeNACA"))
Pkg.add(PackageSpec(path="/Users/markanderson/Box/FLOW-MCA/FLOW-Code/Repositories/personal-projects/NonlinearLiftingLine"))
Pkg.add(PackageSpec(path="/Users/markanderson/Box/FLOW-MCA/FLOW-Code/Repositories/personal-projects/makeAirfoil"))
import NonlinearLiftingLine.NLL
import makeAirfoil.naca

VLM_path = "/Users/markanderson/Box/FLOW-MCA/FLOW-Code/Repositories/personal-projects/VortexLatticeMethod/src/"

include(string(VLM_path,"generatePanels.jl"))

# Straight Wing Geometry
firstCoordinate  = [0.000, 0.000, 0.000];
secondCoordinate = [0.000, 0.500, 0.000];
thirdCoordinate  = [1/6, 0.500, 0.000];
fourthCoordinate = [1/6, 0.000, 0.000];

numPanelsSpan = 100
numPanelsChord = 1 # Must be equal to 1
wingGeometry = generatePanels(firstCoordinate, secondCoordinate, thirdCoordinate, fourthCoordinate, numPanelsSpan, numPanelsChord)
airfoil = naca(0,0,10,0.05)
angleOfAttack = 4.2*pi/180
sideslipAngle = 0
freestream = ones(length(wingGeometry[:,1])).*50

CL, CDi = NLL(wingGeometry,airfoil,angleOfAttack,sideslipAngle,freestream)

println("Nonlinear Lifting Line Theory Results:")
println("CL  = ",CL)
println("CDi = ",CDi)