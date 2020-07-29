using Revise
import CSV
println("3")
# Pkg.instantiate()

import AirfoilPrep
ap = AirfoilPrep
println("8")
import AirfoilDatabase
adb = AirfoilDatabase

import JuliaDB
jdb = JuliaDB

# Database that will be used
database_path = "/Users/markanderson/Desktop/data/practice_database/"

adb.new_database(database_path)
CSV.read(joinpath(database_path, "index.csv"))

#---Creating the data

airfoil_file = "/Users/markanderson/Desktop/data/naca0012.csv"  # Airfoil contour to read
Re = 10e6                           # Reynolds number
Ma = 0.22                           # Mach number
alphas = [i for i in -30:1.0:30]    # Range of AOA to sweep

# Read airfoil contour
x, y = ap.readcontour(airfoil_file; header_len=1, delim=",")

# Run XFOIL generating the polar
polar = ap.runXFOIL(x, y, Re; alphas=alphas,
                              verbose=true, Mach=Ma,
                              iter=100);

# Plot polar
# ap.plot(polar; geometry=true, cdpolar=false);

# Adding polar to database
adb.new_entry(polar; database_path=database_path, airfoilname="NACA 0012")

CSV.read(joinpath(database_path, "index.csv"))

CSV.read(joinpath(database_path, "Cl", "NACA0012-Cl-re10000000-ma0p22-ncrit9p0-0.csv"))

doSweep = false
if doSweep
    alphas = [i for i in -30:1.0:30]    # (deg) angle of attacks to sweep
    Res = [8e4, 1.5e6, 3e6, 5e6, 8e6, 1e7, 1.5e7] # Chord-based Reynolds numbers to sweep
    Mas = [0, 0.15, 0.28]               # Mach numbers to sweep
    ncrits = [14, 9, 4]                 # Ncrit to sweep
    
    # Read airfoil conotur
    x, y = ap.readcontour(airfoil_file; header_len=1, delim=",")
    
    # -------------- RUN SWEEP ------------------------------------------------------
    for ncrit in ncrits
        for Ma in Mas
            for Re in Res
                println("Sweep at Re=$Re\tMa=$Ma\tncrit=$ncrit...")
                
                # Run XFOIL (create polar)
                polar = ap.runXFOIL(x, y, Re; alphas=alphas,
                                              Mach=Ma,
                                              ncrit=ncrit,
                                              verbose=true, 
                                              iter=100)
                # Add polar to the database
                adb.new_entry(polar; database_path=database_path, airfoilname="NACA 0012", warn=false)
                
            end
        end
    end
end

#---RBF Function---#

# Read database
db  = jdb.loadtable(joinpath(database_path, "index.csv"))
Pkg.add("Statistics")
Pkg.add("IterativeSolvers")
import Statistics
import IterativeSolvers

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
    #rbf_scaled, A = generate_RBF(Xp_scaled, val_scaled; zeta=zeta_gaus, sigmas=5.0)
    #rbf(X) = valmin + (valmax - valmin)*rbf_scaled( (X.-Xmin)./ X_scaling )
    
    
    Xps[file] = Xp
    vals[file] = val
    rbfs[file] = rbf
end

display(Xps)
display(vals)
display(rbfs)

xi = 10
testfile = :clfile

println("True value:\t\t$(vals[testfile][xi])")
println("RBF interp value:\t$(rbfs[testfile](Xps[testfile][xi]))")