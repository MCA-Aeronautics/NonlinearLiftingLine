# Purpose: To take airfoil coordinates and use XFOIL to generate coefficient data for it

# include("Desktop/FLOW Lab/Strip Theory/calculateCoefficients.jl")

Pkg.add(PackageSpec(path="/Users/markanderson/Box/FLOW-MCA/FLOW-Code/Repositories/personal-projects/makeAirfoil"))
import makeAirfoil.tabulateData

function calculateCoefficients(chordwiseCoordinates,thicknessCoordinates, angle, reynoldsNumber)

    angle = angle*180/pi # Converting to degrees

    airfoilName = "NACA0010"

    roundedReynoldsNumber = round(reynoldsNumber/100000) * 100000

    filename = string("airfoil-data/",airfoilName,"_",roundedReynoldsNumber,".csv")

    if isfile(filename)
        println(filename," exists!")
        # Then do interpolation
    else
        println(filename," does not exist...")
        println("Let's create that file real fast")
        airfoil = zeros(length(chordwiseCoordinates),2)
        airfoil[:,1] = chordwiseCoordinates
        airfoil[:,2] = thicknessCoordinates
        angleRange = -20:2:30
        path = "airfoil-data"
        airfoilName = "NACA0010"
        makeAirfoil.tabulateData(airfoil,angleRange,path,airfoilName,roundedReynoldsNumber)
    end

    # The old way: each airfoil was calculated on each iteration
    Xfoil.setCoordinates(chordwiseCoordinates,thicknessCoordinates)
    cl, cdd, cdp, cm, converged = Xfoil.solveAlpha(angle,reynoldsNumber)
    #println("Angle = ",angle,", Reynolds Number = ",reynoldsNumber," -> cl = ",cl)
    # return cl[1], cdd[1], cdp[1], cm[1], converged[1]

end