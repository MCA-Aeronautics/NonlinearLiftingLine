# Purpose: To take airfoil coordinates and use XFOIL to generate coefficient data for it

using Pkg
Pkg.add("FLOWMath")
Pkg.add("CSV")
using FLOWMath
using CSV

# Including the tabulateData function
import makeAirfoil.tabulateData

function calculateCoefficients(chordwiseCoordinates,thicknessCoordinates, angle, reynoldsNumber)

    angle = angle*180/pi # Converting to degrees

    airfoilName = "NACA0010"

    filename = string("airfoil-data/",airfoilName,".csv")

    if isfile(filename) == false
        println(filename," does not exist...")
        println("Let's create that file real fast")
        airfoil = zeros(length(chordwiseCoordinates),2)
        airfoil[:,1] = chordwiseCoordinates
        airfoil[:,2] = thicknessCoordinates
        angleRange = -20:1:30
        reynoldsNumberRange = 1e5:1e5:5e6
        path = "airfoil-data"
        tabulateData(airfoil,angleRange,reynoldsNumberRange,path,airfoilName)
    end

    # Do the 2D interpolation
    data = CSV.read(filename)
    data_angles = convert(Array,data[2:end,1])
    data_reynolds = convert(Array,data[1,2:end])
    data_cl = convert(Array,data[2:end,2:end])
    cl = FLOWMath.interp2d(akima,data_angles,data_reynolds,data_cl,angle:angle,reynoldsNumber:reynoldsNumber)

    return cl[1]

end