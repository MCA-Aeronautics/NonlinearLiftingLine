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

    # Eduardo's stuff
    import AirfoilPrep
    ap = AirfoilPrep

    import AirfoilDatabase
    adb = AirfoilDatabase

    import JuliaDB
    jdb = JuliaDB

    Pkg.add("Statistics")
    Pkg.add("IterativeSolvers")
    import Statistics
    import IterativeSolvers


    function NLL(panels,
                 airfoil,
                 airfoilName,
                 freestream,
                 density = 1.225)

        # Defining some variables
        maxIter = 50
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

        # Some of Eduardo's code
        database_path = "/Users/markanderson/Box/FLOW-MCA/Code/NonlinearLiftingLine/airfoil-data/eduardo";
        db  = jdb.loadtable(joinpath(database_path, "index.csv"))

        # Singular kernel
        zeta_sing(r) = r==0 ? 1 : 0

        # erf Gaussian kernel
        zeta_gauserf(r) = 1/(2*pi)^(3/2) * exp(-r^2/2)

        # Gaussian kernel
        zeta_gaus(r) = 3/(4*pi)*exp(-r^3)

        # Winckelmans algebraic kernel
        zeta_wnklmns(r) = 7.5 / (4*pi) / (r^2 + 1)^3.5

        #--MAGIC CODE--#

        # Generate RBF interpolation function
        rbf_axes = [:alpha, :re, :ma, :ncrit]    # Axes of the RBF


        # Read all files in database
        Xps = Dict{Symbol, Array{Array{Float64, 1}, 1}}()
        vals = Dict{Symbol, Array{Float64, 1}}()
        rbfs = Dict{Symbol, Function}()

        # This section just extracts all of the data from the files and feeds it into the generate_RBF() function

        for file in [:clfile, :cdfile, :cmfile]  # Iterate over each file
            
            Xp = Array{Float64, 1}[]
            val = Float64[]
            
            # Convert file symbol to column header
            colheader = Symbol(adb.FIELD2_HEADER[file])
            
            Xmin = [Inf for xi in 1:length(rbf_axes)]
            Xmax = [-Inf for xi in 1:length(rbf_axes)]
            valmin, valmax = Inf, -Inf
            
            for (rowi, row) in enumerate(db) # Iterate over each row

                filename = row[colheader]   # File to read

                # Read data in file
                data = CSV.read(joinpath(database_path, adb.DIRS[file], filename))
                
                for drow in eachrow(data) # Iterate over rows in the data
                    this_Xp = Float64[ax != :alpha ? row[Symbol(adb.FIELD2_HEADER[ax])] : drow[1] for ax in rbf_axes]
                    this_val = drow[2]
                    push!(Xp, this_Xp)
                    push!(val, this_val)
                    
                    for xi in 1:length(rbf_axes)
                        this_Xp[xi] < Xmin[xi] ? Xmin[xi]=this_Xp[xi] : 
                        this_Xp[xi] > Xmax[xi] ? Xmax[xi]=this_Xp[xi] :
                                                nothing
                    end
                    this_val < valmin ? valmin = this_val : this_val > valmax ? valmax = this_val : nothing
                end
                
            end
            
            
            println("Generating $file RBF function with $(length(Xp)) data points...")
            
            # Scale each variable in the range 0 to 1
            X_scaling = [x == 0 ? 1 : x for x in Xmax .- Xmin]
            Xp_scaled = [(X .- Xmin) ./ X_scaling for X in Xp]
            val_scaled = (val .- valmin) ./ (valmax - valmin)
            
            # Generate RBF interpolation function
            rbf, A = generate_RBF(Xp, val; zeta=zeta_gaus, sigmas=0.1)
            
            # Generate scaled RBF interpolation function
        #     rbf_scaled, A = generate_RBF(Xp_scaled, val_scaled; zeta=zeta_gaus, sigmas=1.50)
            rbf_scaled, A = generate_RBF(Xp_scaled, val_scaled; zeta=zeta_gaus, sigmas=5.0)
            rbf(X) = valmin + (valmax - valmin)*rbf_scaled( (X.-Xmin)./ X_scaling )
            
            
            Xps[file] = Xp
            vals[file] = val
            rbfs[file] = rbf
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
                info = [effectiveAOA[j]*180/pi, localReynoldsNumber/chord[j], 0, 9.0]
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
                ylim(minimum(cl_VLM)*0.5,maximum(cl_VLM)*1.25)
                draw()
                println("Iteration: ",i,"\t CL: ",round(CL,digits=6),"\t RMS Difference: ",round(cl_difference_rms,digits=6))

            end

            # Check for convergence
            if cl_difference_rms <= 1*10^-5 # see if the rms difference is small enough to be considered converged           
                break;
            end

            cl_old = cl; # For the next iteration

        end

        # Getting the results from the nonlinear solver
        #CL, cL, _ = calculateLift(density,freestream,panels,GammaValues); # Lift
        CDi_near, _, _ = calculateInducedDrag(density,freestream,panels,GammaValues,cl); # Near-field induced drag

        CL_VLM, cl_VLM, cLSpanLocations_VLM = calculateLift(density,freestream,panels,GammaValues_VLM); # Lift

        CL = trapz(cLSpanLocations_VLM,cl)

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

    # Returns a radial basis function interpolation of a field
    # with values `val` at positions `Xp`. `zeta` is the chosen
    # basis function
    function generate_RBF(Xp, val; zeta=zeta_gaus, sigmas=0.1)
        
        # ERROR CASES
        if size(Xp,1)!=size(val,1)
            error("size(Xp,1)!=size(val,1)")
        end
        
        Np = size(Xp, 1)                     # Number of data points
        
        if size(sigmas)==()
            sgms = sigmas*ones(Np)           # Spreading of every basis function
        else
            sgms = sigmas
        end
            
        # j-th scaled basis evaluated at X
        zetasgm(j, X) = zeta(Statistics.norm(X-Xp[j])/sgms[j])/sgms[j]^3
        
        # Matrix with basis functions evaluated at every point
        # Z[i,j] corresponds to the j-th basis evaluated at i-th point
        Z = [zetasgm(j, Xi) for Xi in Xp, j in 1:Np]
        
        # Solves for the alpha coefficients of every basis
        # A = Z\val
        # A = LinearAlgebra.pinv(Z)*val
        A = IterativeSolvers.gmres(Z, val)
        
        # Generates RBF interpolation function
        rbf(X) = sum([A[j]*zetasgm(j, X) for j in 1:Np])
        
        return rbf, A
    end

end