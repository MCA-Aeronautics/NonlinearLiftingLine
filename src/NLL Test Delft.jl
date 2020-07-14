# This script was used to produce results for the Delft Comparison


using Pkg
Pkg.add("PyPlot")
Pkg.add("CSV")
Pkg.add("FLOWMath")
using PyPlot
using CSV
using FLOWMath

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

numPanelsSpan = 50
numPanelsChord = 1
insidePortion = generatePanels(firstCoordinate,secondCoordinate,thirdCoordinate,fourthCoordinate,numPanelsSpan,numPanelsChord)

# Creating the outer portion of the wing
# Original
firstCoordinate = [0,0.49,0]
secondCoordinate = [0,1.29,0]
thirdCoordinate = [0.161,1.29,0]
fourthCoordinate = [0.279,0.49,0]

# Experimental
# firstCoordinate = [0,0.49,0]
# secondCoordinate = [0.161/4,1.29,0]
# thirdCoordinate = [0.161,1.29,0]
# fourthCoordinate = [0.279,0.49,0]

numPanelsSpan = 50

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

CL, CDi, cl, spanLocations = NLL(sortedPanels,airfoil,"NACA642-015A",freestream)

# Get the Delft VLM and experimental data
data = CSV.read("/Users/markanderson/Box/FLOW-MCA/Lab Notebooks/2020-07/Delft VLM Data.csv")
spancoords = convert(Array,data[1:end,1])
normalizedLift = convert(Array,data[1:end,2])

# Calculating the chord
chord = zeros(numPanels,1)
for i = 1:numPanels
    chord_lhs = sortedPanels[i,10] - sortedPanels[i,1]
    chord_rhs = sortedPanels[i,7] - sortedPanels[i,4]
    chord[i] = (chord_lhs + chord_rhs) / 2
end

figure()
plot(spanLocations./maximum(spanLocations),cl.*chord./0.24,label="Strip Theory", color = "orange", linestyle = "--", linewidth = 3)
scatter(spancoords,normalizedLift,label="Delft VLM", color = "black")
xlim(0,1)
title("Normalized CL Comparison")
legend()
xlabel("Span Location (2y/b)")
ylabel("cl * chord / 0.24")

# Creating an unnormalized lift coefficient array for the Delft data
unnormalizedLift = zeros(length(normalizedLift),1)
for i = 1:length(unnormalizedLift)
    unnormalizedLift[i] = normalizedLift[i] * 0.24 / linear(spanLocations./maximum(spanLocations),chord,spancoords[i])
end

figure()
plot(spanLocations./maximum(spanLocations),cl,label="Strip Theory", color = "orange", linestyle = "--", linewidth = 3)
scatter(spancoords,unnormalizedLift,label="Delft VLM",color = "black")
xlim(0,1)
title("CL Comparison")
legend()
xlabel("Span Location (2y/b)")
ylabel("cl")