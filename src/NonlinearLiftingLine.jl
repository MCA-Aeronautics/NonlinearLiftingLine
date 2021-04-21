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
    # Pkg.add("PyPlot")
    # Pkg.add("LinearAlgebra")
    # Pkg.add("FLOWMath")
    # Pkg.add("CSV")
    # Pkg.develop(PackageSpec(url="https://github.com/byuflowlab/Xfoil.jl"))
    # Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/VLMMCA"))
    # Pkg.develop(PackageSpec(url="https://github.com/Mark-C-Anderson/makeAirfoil"))
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

    # Eduardo's stuff
    import AirfoilPrep
    ap = AirfoilPrep

    import AirfoilDatabase
    adb = AirfoilDatabase

    import JuliaDB
    jdb = JuliaDB

    #Pkg.add("Statistics")
    #Pkg.add("IterativeSolvers")
    import Statistics
    import IterativeSolvers


    function NLL(panels,
                 airfoil,
                 airfoilName,
                 freestream,
                 rbfs,
                 density = 1.225)

        # Defining some variables
        maxIter = 1000
        showConvergence = false
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
        #GammaValues = ones(numPanels,1)*0.01

        # Initializing the coefficient arrays
        cl = zeros(numPanels) # section lift coefficient
        cdd = zeros(numPanels) # section induced drag coefficient
        cdp = zeros(numPanels) # section pressure drag coefficient
        cm = zeros(numPanels)  # section moment coefficient
        converged = zeros(numPanels) # whether the iteration converged

        # Calculating the chord
        chord = zeros(numPanels,1)
        for i = 1:numPanels
            chord_lhs = panels[i,10] - panels[i,1]
            chord_rhs = panels[i,7] - panels[i,4]
            chord[i] = (chord_lhs + chord_rhs) / 2
        end

        # Now we will create a loop that will converge to a circulation value that will be available
        # to be used to calculate our aerodyanmic coefficients and forces
        for i = 1:maxIter # Maximum number of iterations to go through before giving up (to avoid an infinite loop)

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

            # if i <= 5
            #     figure(3)
            #     plot(spanLocations,GammaValues)
            #     title("Circulation Values")
            # end

            # if i <= 5
            #     figure(4)
            #     plot(spanLocations,inducedVelocity)
            #     title("Induced Velocity")
            # end

            # if i <= 5
            #     figure(5)
            #     plot(spanLocations,effectiveAOA)
            #     title("Effective AOA")
            # end

            # find the total circulation
            for j = 1:(numPanels)

                # accounting for the difference velocity at each airfoil
                # FIXME: Why only in the x-direction?
                #localVelocity = sqrt(freestream[j,1]^2 + inducedVelocity[j]^2)
                localVelocity = sqrt(freestream[j,1]^2 + (freestream[j,3] + inducedVelocity[j])^2)

                localReynoldsNumber = localVelocity * chord[j] / nu

                # Calculate the lift and drag coefficients for that angle for each airfoil
                #cl[j] = calculateCoefficients(airfoil[:,1],airfoil[:,2], effectiveAOA[j], localReynoldsNumber/chord[j], airfoilName);
                info = [effectiveAOA[j]*180/pi, localReynoldsNumber, 0, 9.0]
                cl[j] = rbfs[:clfile](info)

                # println(info)
                # println(cl[j])
                
                #println("Iteration ",i," Panel ",j,". cl = ",cl[j])

                # Use the coefficients to calculate the circulation about each airfoil.
                # When you define the local lift coefficient as cl = L / (q * S) and replace L with Kutta-Joukowski,
                # you can simplify to the following expression. Remember that deltaX = chord.
                circulation = 0.5 * norm(freestream[1,:]) * cl[j] * chord[j]

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
                # plot(spanLocations,cl_VLM.*chord./0.24,label = "Linear Result",color = "green")
                # plot(spanLocations,cl.*chord./0.24,label = "Nonlinear Result",linestyle = "--",color = "orange",linewidth = 3)
                xlabel("Spanwise Location (y/b)")
                ylabel("Sectional Lift Coefficient")
                title(string("Iteration ",i))
                grid("on")
                legend()
                #xlim(-1.3,1.3)
                ylim(minimum(cl)*0.5,maximum(cl)*1.25)
                draw()
                println("Iteration: ",i,"\t CL: ",round(CL,digits=6),"\t RMS Difference: ",round(cl_difference_rms,digits=6))

            end

            # Check for convergence
            if cl_difference_rms <= 1*10^-4 # see if the rms difference is small enough to be considered converged           
                break;
            end

            cl_old = cl; # For the next iteration

        end

        # Getting the results from the nonlinear solver
        CL, cL, _ = calculateLift(density,freestream,panels,GammaValues); # Lift
        CDi_near, _, _ = calculateInducedDrag(density,freestream,panels,GammaValues,cl); # Near-field induced drag

        CL_VLM, cl_VLM, cLSpanLocations_VLM = calculateLift(density,freestream,panels,GammaValues_VLM); # Lift

        #CL = trapz(cLSpanLocations_VLM,cl)

        return CL, CDi_near, cl, cLSpanLocations_VLM

    end

    function calculateCoefficients(database_path)

        angle = angle*180/pi # Converting to degrees
    
        # airfoilName = "NACA0010"
    
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