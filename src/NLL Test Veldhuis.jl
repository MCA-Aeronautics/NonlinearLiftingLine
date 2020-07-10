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
secondCoordinate = [0,0.640,0]
thirdCoordinate = [0.240,0.640,0]
fourthCoordinate = [0.240,0,0]

# Experimental
# firstCoordinate = [0,0,0]
# secondCoordinate = [0,0.5,0]
# thirdCoordinate = [0.2,0.5,0]
# fourthCoordinate = [0.2,0,0]

numPanelsSpan = 100
numPanelsChord = 1
panels = generatePanels(firstCoordinate,secondCoordinate,thirdCoordinate,fourthCoordinate,numPanelsSpan,numPanelsChord)

sortedPanels = sortslices(panels,dims=1,by=x->x[2],rev=false) # sorted by the firstCoordinate y-value

numPanels = length(panels[:,1])

freestream = zeros(numPanels,3)

Vinf = 50

alpha = 4*pi/180

for i = 1:numPanels
    freestream[i,:] = [Vinf*cos(alpha),0,Vinf*sin(alpha)]
end

#------------
# Getting the wake data ready
#------------

wakeData = CSV.read("Propeller Data/E212 Wake.csv")

propDiameter = 0.236 # meters

propLocation = 0.300 # meters

radialPosition = wakeData[:,1]
axialVelocity = wakeData[:,2]
swirlVelocity = wakeData[:,3]

# Determine whether each panel is within the wake. If it is, then add the axial and swirl components to the freestream
for i = 1:numPanels

    # determine distance to center of Propeller
    distanceFromProp = abs(panels[i,2]) - propLocation

    if abs(distanceFromProp) <= propDiameter/2 # See if it's within the wake
        axial = akima(radialPosition,axialVelocity,abs(distanceFromProp))
        swirl = akima(radialPosition,swirlVelocity,abs(distanceFromProp))

        # If you are to the outside of the propeller, reverse the swirl because the propeller is
        # swinging down
        if distanceFromProp > 0
            swirl = -swirl
        end

        freestream[i,:] = freestream[i,:] + [axial, 0, swirl]

    end

end



# modelProp = 0:(2*pi/30):(2*pi)
# modelWash = similar(modelProp)
# for i = 1:length(modelProp)

#     modelWash[i] = 5 * sin(modelProp[i])

# end

# freestream[35:65,:] = freestream[35:65,:] .+ modelWash

# freestream[135:165,:] = freestream[135:165,:] .+ reverse(modelWash)


#-----------
# On with the calculations
#-----------

pygui(true)

CL, CDi, cl, spanLocations = NLL(sortedPanels,airfoil,"NACA642-015A",freestream)

# Get the Delft VLM and experimental data
data = CSV.read("/Users/markanderson/Box/FLOW-MCA/Lab Notebooks/2020-07/Veldhuis Propeller Data.csv")
spancoords = convert(Array,data[1:end,1])
unnormalizedLift = convert(Array,data[1:end,2])

figure()
plot(spanLocations./maximum(spanLocations),cl,label="Strip Theory", color = "orange", linestyle = "--", linewidth = 3)
plot(spancoords./maximum(spanLocations),unnormalizedLift,label="Veldhius Data", color = "black",marker = "o")
xlim(0,1)
title("Veldhuis CL Comparison")
legend()
xlabel("Span Location (2y/b)")
ylabel("cl")
