# ADRE-Reactive-Transport-Model
This reactive transport model contains Julia code for the numerical simulation and optimization of a 2 phase mass transport and reaction process. The model simulates the addition of a treatment method (H2O2) in a bio clogged 1D porous media. The first phase simulates biofilm decay, while the second phase simulates growth.

The model takes various inputs including numerical discretization, mass transfer parameters, reaction kinetic parameters, and porosity-permeability parameters.

The model begins with numerical simutation of biofilm growth and decay across a two phase period. Using initial and intermediary conditions, the ADRE is solved using a numerical differential equation solver.

The ADRE solves simoultaneosly the transport of H2O2 and solute, and the decay and growth of biofilm.
Reaction kinetics for the growth of bifilm are simulated using a monod model.
Decay kinetics are simulated as first order interactions with H2O2 free radicals.

The primary output of the ADRE simulation is biofilm volume fraction.
Whithin the numerical simulation an exponential model is applied to convert biofilm volume fraction to biofilm porosity, and a second exponential model is used to convert biofilm volume fraction to column permeability.
Alpha and Beta are target variables in the exponential functions that are fitted

Numerical Simutlation and conversion to permeability are contained in a block called sim.function().
This sim.function() is then utilized with a sum of least squares error Nelder Mead optimization algorithm with Alpha and Beta to fit the parameters to experimental data.

A plot is produced visualizing experimental data and subsequent curve fit. 
