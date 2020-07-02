# Purpose: Model a wing by solving for aerodynamic properties at each airfoil segment
# Inputs:
#
# Outputs:
#
# TIPS:
# The freestream should be larger than 1m/s because XFOIL will report separation for these values. Your reynolds number should be at least 100,000

module NonlinearLiftingLine

    # Adding all the additional packages that are used in this code
    using Pkg
    Pkg.add("PyPlot")
    Pkg.add("LinearAlgebra")
    Pkg.add("FLOWMath")
    Pkg.add("CSV")
    Pkg.develop(PackageSpec(url="https://github.com/byuflowlab/Xfoil.jl"))
    Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
    Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
    using PyPlot
    using LinearAlgebra
    using FLOWMath
    using CSV
    import Xfoil
    import VLMMCA.VLM
    import VLMMCA.calculateLift
    import VLMMCA.calculateInducedDrag
    import VLMMCA.calculateInducedVelocity
    #import VLMMCA.
    import makeAirfoil.tabulateData
    #import makeAirfoil.


    function NLL(panels,
                 airfoil,
                 freestream,
                 density = 1.225)

        # Defining some variables
        maxIter = 500
        showConvergence = true
        numPanels = length(panels[:,1])
        nu = 1.48e-5;
        cl_old = zeros(numPanels);

        # Calculating the angles of attack across the span
        anglesOfAttack = zeros(numPanels)
        for i = 1:numPanels
            anglesOfAttack[i] = atan(freestream[i,3] ./ freestream[i,1])
        end

        CL_VLM, CDi_near_VLM, cl_VLM, cd_near_VLM, spanLocations, GammaValues_VLM = VLM(panels,freestream,density)

        # Initialize the gamma values
        GammaValues = GammaValues_VLM*1.0

        # Initializing the coefficient arrays
        cl = zeros(numPanels) # section lift coefficient
        cdd = zeros(numPanels) # section induced drag coefficient
        cdp = zeros(numPanels) # section pressure drag coefficient
        cm = zeros(numPanels)  # section moment coefficient
        converged = zeros(numPanels) # whether the iteration converged

        # Now we will create a loop that will converge to a circulation value that will be available
        # to be used to calculate our aerodyanmic coefficients and forces
        for i = 1:maxIter + 1 # Maximum number of iterations to go through before giving up (to avoid an infinite loop)

            # Initializing the coefficient arrays on each iteration
            cl = zeros(numPanels) # section lift coefficient
            cdd = zeros(numPanels) # section induced drag coefficient
            cdp = zeros(numPanels) # section pressure drag coefficient
            cm = zeros(numPanels)  # section moment coefficient
            converged = zeros(numPanels) # whether the iteration converged
            totalCirculation = 0.0 # initialized to zero
            oldGammaValues = GammaValues

            # Calculate the induced velocity at the front of each horseshoe vortex (1/4 chord)
            inducedVelocity = calculateInducedVelocity(panels,GammaValues,"quarter chord") 

            # Calculate the effective angle of attack for each airfoil
            effectiveAOA = calculateEffectiveAlpha(freestream,inducedVelocity,anglesOfAttack) # multiplied by cosine of the angle of attack so that it becomes perpendicular to the freestream
            
            # find the total circulation
            for j = 1:(numPanels)

                # accounting for the difference velocity at each airfoil
                localVelocity = sqrt(dot(freestream[j,:],[1,0,0])^2 + inducedVelocity[j]^2)
                chord = panels[j,10] - panels[j,1]
                localReynoldsNumber = localVelocity * chord / nu

                # Calculate the lift and drag coefficients for that angle for each airfoil
                cl[j] = calculateCoefficients(airfoil[:,1],airfoil[:,2], effectiveAOA[j], localReynoldsNumber/chord);
                
                # Use the coefficients to calculate the circulation about each airfoil. Not sure which one to use
                circulation = 0.5 * freestream[j] * cl[j] * chord

                # update the GammaValues array
                omega = 0.99
                GammaValues[j] = omega*oldGammaValues[j] + (1-omega)*circulation[1]

            end

            # Defining the rms difference between the previous cl and current cl distributions
            cl_difference = cl .- cl_old
            cl_difference_rms = norm(cl_difference)/sqrt(length(cl_difference))

            # Indicate how far along it is and plot the distribution
            if showConvergence

                CL,_,_ = calculateLift(density,freestream,panels,GammaValues); # Lift

                figure(1)
                clf()
                plot(spanLocations,cl_VLM,label = "Linear Result",color = "green")
                plot(spanLocations,cl,label = "Nonlinear Result",linestyle = "--",color = "orange",linewidth = 3)
                xlabel("Spanwise Location (y/b)")
                ylabel("Sectional Lift Coefficient")
                title(string("Iteration ",i))
                grid("on")
                legend()
                xlim(-0.5,0.5)
                ylim(minimum(cl_VLM)*0.5,maximum(cl_VLM)*1.25)
                draw()
                println("Iteration: ",i,"\t CL: ",round(CL,digits=6),"\t RMS Difference: ",round(cl_difference_rms,digits=6))

            end

            # Check for convergence
            if cl_difference_rms <= 1*10^-6 # see if the rms difference is small enough to be considered converged           
                break;
            end

            cl_old = cl; # For the next iteration

        end

        # Getting the results from the nonlinear solver
        CL, cl, _ = calculateLift(density,freestream,panels,GammaValues); # Lift
        CDi_near, _, _ = calculateInducedDrag(density,freestream,panels,GammaValues,cl); # Near-field induced drag

        CL, cl_VLM, cLSpanLocations_VLM = calculateLift(density,freestream,panels,GammaValues_VLM); # Lift

        return CL, CDi_near, cl

    end

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
    
    function calculateEffectiveAlpha(freestream,inducedVelocity,anglesOfAttack)

        inducedAlpha = zeros(length(freestream[:,1]))
        effectiveAlphas = similar(inducedAlpha)
        for i = 1:length(freestream[:,1])
            
            inducedAlpha[i] = atan(inducedVelocity[i],freestream[i,1]) # Just x-component of freestream
    
            effectiveAlphas[i] = anglesOfAttack[i] .+ inducedAlpha[i]
    
        end
    
        return effectiveAlphas
    
    end

end