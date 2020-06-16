revise()
Pkg.develop(PackageSpec(path=pwd()))
import NonlinearLiftingLine.NLL

include("../VortexLatticeMethod/src/generatePanels.jl")

include("../makeAirfoil/src/makeAirfoil.jl")
import makeAirfoil.naca

# Straight Wing Geometry
firstCoordinate  = [0.000, 0.000, 0.000];
secondCoordinate = [0.000, 0.500, 0.000];
thirdCoordinate  = [1/6, 0.500, 0.000];
fourthCoordinate = [1/6, 0.000, 0.000];

numPanelsSpan = 100
numPanelsChord = 1 # Must be equal to 1
numPanels = numPanelsSpan*numPanelsChord*2
wingGeometry = generatePanels(firstCoordinate, secondCoordinate, thirdCoordinate, fourthCoordinate, numPanelsSpan, numPanelsChord)
airfoil = naca(0,0,10,0.05)
angleOfAttack = 4.2*pi/180
sideslipAngle = 0
freestream = zeros(numPanels,3)
for i = 1:numPanels
    freestream[i,:] = 50 .* [cos(angleOfAttack)*cos(sideslipAngle),-sin(sideslipAngle),sin(angleOfAttack)*cos(sideslipAngle)]
end

CL, CDi = NLL(wingGeometry,airfoil,freestream)

println("Nonlinear Lifting Line Theory Results:")
println("CL  = ",CL)
println("CDi = ",CDi)