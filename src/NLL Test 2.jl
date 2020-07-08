using Pkg
Pkg.add("PyPlot")
Pkg.add("CSV")
using PyPlot
using CSV

revise()

# Getting the VLMMCA package
Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
#Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/NonlinearLiftingLine"))
Pkg.develop(PackageSpec(path=pwd()))

import VLMMCA.generatePanels
import VLMMCA.plotPanels
import VLMMCA.VLM

import makeAirfoil.naca

import NonlinearLiftingLine.NLL

# Creating the airfoil
# airfoil = naca(0,0,20,0.05)

data = CSV.read("/Users/markanderson/Box/FLOW-MCA/Lab Notebooks/2020-07/Delft-airfoil.csv")
xcoords = convert(Array,data[1:end,1])
ycoords = convert(Array,data[1:end,2])
airfoil = cat(xcoords,ycoords,dims=2)

# Creating the inner portion of the wing
# Original
firstCoordinate = [0,0,0]
secondCoordinate = [0,0.49,0]
thirdCoordinate = [0.279,0.49,0]
fourthCoordinate = [0.279,0,0]

# Experimental
# firstCoordinate = [0,0,0]
# secondCoordinate = [0,0.5,0]
# thirdCoordinate = [0.2,0.5,0]
# fourthCoordinate = [0.2,0,0]

numPanelsSpan = 100
numPanelsChord = 1
insidePortion = generatePanels(firstCoordinate,secondCoordinate,thirdCoordinate,fourthCoordinate,numPanelsSpan,numPanelsChord)

# Creating the outer portion of the wing
# Original
firstCoordinate = [0,0.49,0]
secondCoordinate = [0,1.29,0]
thirdCoordinate = [0.161,1.29,0]
fourthCoordinate = [0.279,0.49,0]

# Experimental
# firstCoordinate = [0,0.5,0]
# secondCoordinate = [0.1,1.0,0]
# thirdCoordinate = [0.3,1.0,0]
# fourthCoordinate = [0.2,0.5,0]

numPanelsSpan = 100

outsidePortion = generatePanels(firstCoordinate,secondCoordinate,thirdCoordinate,fourthCoordinate,numPanelsSpan,numPanelsChord)

panels = cat(insidePortion,outsidePortion,dims=1)

sortedPanels = sortslices(panels,dims=1,by=x->x[2],rev=false) # sorted by the firstCoordinate y-value

numPanels = length(panels[:,1])

freestream = zeros(numPanels,3)

Vinf = 19

alpha = 4*pi/180

for i = 1:numPanels
    freestream[i,:] = [Vinf*cos(alpha),0,Vinf*sin(alpha)]
end

pygui(true)

CL, CDi, cl = NLL(sortedPanels,airfoil,"NACA642-015A",freestream)