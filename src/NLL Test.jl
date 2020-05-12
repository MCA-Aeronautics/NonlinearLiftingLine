# To include this file:
# include("Box/FLOW Lab/Modules MCA/NonlinearLiftingLine/src/NLL Test.jl")

using Pkg; Pkg.add("Revise"); using Revise;
Pkg.develop(PackageSpec(path = "/Users/markanderson/Box/FLOW Lab/Modules MCA/makeNACA"))
Pkg.develop(PackageSpec(path = "/Users/markanderson/Box/FLOW Lab/Modules MCA/NonlinearLiftingLine/"))
import makeNACA.naca
import NonlinearLiftingLine.NLL

VLM_path = "/Users/markanderson/Box/FLOW Lab/Modules MCA/VortexLatticeMethod/src/"

include(string(VLM_path,"generatePanels.jl"))

# Straight Wing Geometry
firstCoordinate  = [0.000, 0.000, 0.000];
secondCoordinate = [0.000, 0.500, 0.000];
thirdCoordinate  = [1/6, 0.500, 0.000];
fourthCoordinate = [1/6, 0.000, 0.000];

numPanelsSpan = 20
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