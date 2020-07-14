revise()
Pkg.develop(PackageSpec(path=pwd()))
Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
import NonlinearLiftingLine.NLL
import VLMMCA.generatePanels
import makeAirfoil.naca

# Straight Wing Geometry
firstCoordinate  = [0.000, 0.000, 0.000];
secondCoordinate = [0.000, 0.500, 0.000];
thirdCoordinate  = [1/6, 0.500, 0.000];
fourthCoordinate = [1/6, 0.000, 0.000];

# # Swept Wing Geometry
# firstCoordinate  = [0.000, 0.000, 0.000];
# secondCoordinate = [1.0, 1.0, 0.000];
# thirdCoordinate  = [1.2, 1.0, 0.000];
# fourthCoordinate = [0.2, 0.000, 0.000];

numPanelsSpan = 50
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

# modelProp = 0:(2*pi/30):(2*pi)
# modelWash = similar(modelProp)
# for i = 1:length(modelProp)

#     modelWash[i] = 5 * sin(modelProp[i])

# end

# freestream[35:65,:] = freestream[35:65,:] .+ modelWash

# freestream[135:165,:] = freestream[135:165,:] .+ reverse(modelWash)

CL, CDi, _ = NLL(wingGeometry,airfoil,freestream)

println("Nonlinear Lifting Line Theory Results:")
println("CL  = ",CL)
println("CDi = ",CDi)