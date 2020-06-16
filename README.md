# NonlinearLiftingLine
A strip theory code that I wrote as a member of the BYU FLOW Lab

# Getting started
<pre><code> git clone --recursive https://github.com/Mark-C-Anderson/NonlinearLiftingLine.git
</code></pre>

## Instructions for running the code
Make sure that you are in the NonlinearLiftingLine directory, but above the src/ folder. Open a Julia REPL and type:

<pre><code> include("src/NLL Test.jl") </code></pre>

# How this code works
This is a "strip theory" code, but the more proper name would be "Nonlinear Lifting Line". It works by:

1. Taking in a wing geometry with an arbitrary number of spanwise panels and a single chordwise panel. A single airfoil is used across the wing.
2. If data are not yet generated for a given airfoil, those data will be generated. As of this release, the default name is "NACA0010". This may be updated to take user input in the future. Data are stored in the airfoil-data/ directory as CSV files.
3. An initial constant circulation distribution is assummed, and the panels are each treated as horseshoe vortices. The downwash induced at each panel is then calculated using this distribution of horseshoe vortices.
4. The effective angle of attack is calculated at each panel, and this along with the local Reynolds number is used to interpolate on the airfoil data to return an updated lift coefficient value.
5. The updated lift coefficient values are used to calculate a new circulation distribution.
6. This process is iterated until changes in the circulation are small (<1e6).
7. Results are plotted on each iteration and compared with the results from the Linear Lifting Line method (Vortex Lattice Method)
