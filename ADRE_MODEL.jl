### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# ╔═╡ 3c915980-5601-11f0-34f2-37291524e5a4
begin
	import Pkg; Pkg.activate()
	using Controlz, CairoMakie, DataFrames, Optim, MakieThemes, Dates, CSV, DifferentialEquations, Printf, Statistics, LinearAlgebra, MAT, Optim, Interpolations, Latexify


	set_theme!(ggthemr(:fresh))
	update_theme!(fontsize=24, markersize=16, linewidth=3)
end

# ╔═╡ ca47908b-cfef-4ced-ab9f-3e8d89e69970
#1-D Reactive Transport Model with Two Phases :
#(Phase 1) H2O2 Injection (0 ≤ t ≤ 50 [hrs])
#(Phase 2) Nutrient Injection (Treated Wastewater) (50 < t ≤ 120 [hrs])

#Biofilm fraction ϕ_b decays in Phase (A) due to H2O2 flux, then regrows in Phase (B) due to lack of H2O2 flux.

#For this model, we track:
#-C(x,t): H2O2 concentration (mg/l) using ADRE mass transfer
#-ϕ_b(x,t): biofilm fraction [0 ≤ ϕ_b ≤ 1] using monod kinetics
#-K(x,t): permeability using exponential relationships 


#Mads Kent -- July 1, 2025


# ╔═╡ 8c4ee984-9339-4f9e-af82-66ae9c342ce6
begin #1) MODEL PARAMETERS

#Geometry and Discretization
	const L = 0.3 #Column Length [m]
	const N = 60 #Discretization intervals []
	const dx = L / N #[m]
	x = range(0,L,N+1)
	const A_column = 0.002 #[m²]
	
end

# ╔═╡ a417ee46-af67-4d9e-af56-d3bdb621a2c9
begin #Flow & Time
	
	const Q_in = 2.25e-3 / 24.0 #[m³/hr]
	const u_0 = Q_in / A_column #Darcy velocity[m/s] 
	
end

# ╔═╡ 53beaddf-502a-4464-8891-0a4928739e2f
begin #Simulate Experiment Until 120 hrs [0 - 50 is H2O2 injection]
	
	const T_inject_H2O2 = 50.0 #H2O2 injection time [hr]
	const T_total = 120.0 #Total simulation time [hr]
	
end

# ╔═╡ 312fa5d6-4a9d-4154-95a8-66820e555fea
begin #Porosity & Initial Biofilm
	
	const ϕ_0 = 0.349 #from micromodel [] 
	const χ_0 = 0.69 #Biofilm volume fraction [] 
	
end

# ╔═╡ 750a7e52-9e4c-4eb9-9a68-622dca3ae870
begin #Transport Parameters
	
	const D_c = 3.0e-5 #Diffusion Coefficient for H2O2 [m²/hr]
	const D_s = 2.5e-5 #Diffusion Coefficient for substrate [m²/hr]
	
end

# ╔═╡ d883de23-c596-423b-8866-e457a232a8d9
begin #H2O2 Reaction/Decay
	
	const k_ox = 5.0e-4 #H2O2 consumption rate [1/hr]
	const k_death = 2.25e-4 #Oxidative biofilm decay rate
	
end

# ╔═╡ b2525aea-5c06-416a-9efe-dfb34bfbdda9
begin #Monod growth for phase 2
	
	const μ_max = 3.75e-2 #max growth rate [1/hr]
	const K_s = 40.0 #half-sat constant [mg/l]
	const Y = 10.0 #yield coefficient []
	
end

# ╔═╡ f9dadcaa-72e8-4383-afc5-cca26b6a5673
begin #Inlet Boundary Conditions
	
	const C_in_phase1 = 50.0 #[mg/l]
	const S_in_phase1 = 0.0 #[mg/l]
	const C_in_phase2 = 0.0 #[mg/l]
	const S_in_phase2 = 30.0 #[mg/l]
	
end

# ╔═╡ f5ab7022-a472-414c-99c0-9e9f0181382c
begin #Porosity Relationship
	
	α = 1.5 #
	
end

# ╔═╡ 12c3b492-8f63-4b3e-a1d7-d0f346865eea
begin #Permeability Relationship
	
	K_0 = 1.0793e-13 #[m²] Permeability when biofilm is already in the system
	β = 1.2 #[] (1.2)
	B_c = 0.6
	K_bio = 1.05e-12
	
end

# ╔═╡ b11fc159-2eaf-46d4-ac65-c5a40ab2ed6a
#Define general exponential function for use in porosity and permeability models.
	function exponential(A_0, τ, t)
		T = A_0 * exp(-τ * t)
		return T
	end

# ╔═╡ 719d7bbf-d6b2-49e1-b5cf-f9ed59f40953
begin # 2) INITIAL FIELDS

	#H2O2 Array
	C0 = zeros(N+1)

	#Substrate Array
	S0 = zeros(N+1)

	#Biofilm fraction
	χ0 = χ_0 * ones(N+1)

end

# ╔═╡ a8130b51-b286-4768-96ff-949dcccc55df
begin
function reactor!(du, u, p, t)
	C, S, χ = @views u[1:N+1], u[N+2:2(N+1)], u[2(N+1)+1:end]
	du_C, du_S, du_χ = @views du[1:N+1], du[N+2:2(N+1)], du[2(N+1)+1:end]

	#Inlet Boundary Conditions
	C_in = t <= T_inject_H2O2 ? C_in_phase1 : C_in_phase2
    S_in = t <= T_inject_H2O2 ? S_in_phase1 : S_in_phase2
    C[1] = C_in
    S[1] = S_in

	#H2O2 Transport
	@inbounds for i in 2:N
		adv_flux = -(u_0/dx) * (C[i] - C[i-1])
		diff_flux = D_c * (C[i+1] - 2C[i] + C[i-1]) / dx^2
		react_term = k_ox * χ[i] * C[i]
		du_C[i] = (adv_flux + diff_flux - react_term) / ϕ_0
	end

	C[N+1] = C[N]
	du_C[N+1] = 0.0
	
	#Substrate Transport
	@inbounds for i in 2:N
        adv_flux_s = -(u_0/dx) * (S[i] - S[i-1])
        diff_flux_s = D_s * (S[i+1] - 2S[i] + S[i-1]) / dx^2
        growth_term = μ_max * (S[i]/(K_s + S[i])) * χ[i]
        consume_term = Y * growth_term
        du_S[i] = (adv_flux_s + diff_flux_s - consume_term) / ϕ_0
	end
	
	S[N+1] = S[N]
	du_S[N+1] = 0.0
	
	#Biofilm Kinetics
	@inbounds for i in 1:N+1
        growth_i = μ_max * (S[i]/(K_s + S[i])) * χ[i]
        decay_i = k_death * C[i] * χ[i]
        du_χ[i] = growth_i - decay_i
        # Ensure phi_b ∈ [0, 1]
        χ[i] = clamp(χ[i], 0.0, 1.0)
    end	
end
end

# ╔═╡ a87a36eb-fd73-464a-9177-3b6a18fd7759


# ╔═╡ 8a3cd7ea-29fc-46be-b74f-26a282166017
begin # This function is the target function to use in the Nelder Mead Optimization
function sim_function(α, β, K_optim, T_0, T_f)

	#Prepare numerical simulation
	u0 = vcat(C0, S0, χ0)
	tspan = (0.0, T_total)
	prob = ODEProblem(reactor!, u0, tspan)
	
	sol = solve(prob, Tsit5(), saveat=0.1)

	# Define time range   

	t_idx = findall(t -> T_0 ≤ t ≤ T_f, sol.t)
	
	χ_range = sol[2*(N+1)+1:end, t_idx]

	t_range = sol.t[t_idx]
	
	ϕ_field = exponential.(ϕ_0, α, χ_range)# Porosity Model
	
	K_field = exponential.(K_optim, -β, ϕ_field) # Permeability Model
	
	K_avg = vec(mean(K_field, dims=1))
	
	K_avg_norm = K_avg / K_0 #Normalized permeability
	
	return K_avg_norm, sol, t_range
end
end

# ╔═╡ 3f436653-0912-423b-b087-7888ccbef80c
begin # Retrieve experimental data, specfically Matlab data.
	file = matopen("C:/Users/madsc/Downloads/Perez_Kdata_H2O2treatment.mat")
	varnames = keys(file)
	println(varnames)

	t_data = vec(read(file, "tH2O2start"))
	normK_data = vec(read(file, "normK"))
	close(file)	
end

# ╔═╡ 357ea7cb-9e0a-4e6e-98f2-f124fcedf5fb
function loss_coupled(x) # Define the loss function that calculates SSE
    α_coupled, β_coupled, K_optim = x # prepare parameter vector
	
    K_sim, sol, t_range = sim_function(α_coupled, β_coupled, K_optim, 0.0, T_total)
    # Ensure K_sim is a vector with same length as sol.t
    K_sim_vec = vec(K_sim)

    itp = linear_interpolation(t_range, K_sim_vec, extrapolation_bc=Line())
    
    sse = 0.0
    for (t_obs, K_obs) in zip(t_data, normK_data)
        K_pred = itp(t_obs)
        sse += (K_obs - K_pred)^2
    end
    
    return sse
end

# ╔═╡ 469194e5-0e01-439a-9cca-8dcbec018c1f
begin # Initial guesses for parametric optimization. 
	α_coupled_intital = 1.5

	β_coupled_intital = 1.4

	K_optim_coupledinitial = 1.0793e-13
	
	initial_params_coupled = [α_coupled_intital, β_coupled_intital, K_optim_coupledinitial]
	
	lower_coupled = [-Inf, -Inf, -Inf] #Useful to bind values 

	upper_coupled = [Inf, Inf, Inf]
end

# ╔═╡ 2452d35f-f2b7-4202-a989-decbc58b260d
begin 
	#Perform Nelder Mead optimization with SSE loss function
	result_coupled = optimize(loss_coupled, lower_coupled, upper_coupled, initial_params_coupled, Fminbox(NelderMead()), Optim.Options(iterations=1000))
println("Optimized parameters: ", Optim.minimizer(result_coupled))
end

# ╔═╡ 07bd9879-d8d8-4c27-94c4-a6ce4437ba55
begin
	#extract optimization result
	α_coupledoptim, β_coupledoptim, K_optim_coupled_optim = Optim.minimizer(result_coupled)
	
	K_avg_norm_optim_coupled, sol_optim_coupled, t_range_optim_coupled = sim_function(α_coupledoptim, β_coupledoptim, K_optim_coupled_optim, 0.0, T_total)

	Latex_data = [t_range_optim_coupled;K_avg_norm_optim_coupled]
	
	t_range_optim_coupled_pv = t_range_optim_coupled.*( Q_in / (A_column * L * ϕ_0))
	
end

# ╔═╡ d4a9accb-d00d-459d-ac7b-e15ea50dd71d
begin
	# Format figure and plot data
	Permfig = Figure()
	ax = Axis(Permfig[1,1],
	title = "Biofilm Permeability Time Series Data",
	xlabel = "Time [hrs]",
	ylabel = "Permeability of Biofilm [m²] ",
	titlefont = :bold, 
    titlesize = 18,
    xlabelsize = 14,
    ylabelsize = 14,
	limits=(0, 130, 0, 2) ) 

	scatter!(ax, t_data, normK_data, 
	markersize = 10,        
    marker = :rect,         
    strokewidth = 2,        
    strokecolor = :green,   
    color = (:transparent),
	label = "Empirical Data"
	)
	
	lines!(ax, t_range_optim_coupled, K_avg_norm_optim_coupled,
		label = "One Alpha and Beta")
	
	axislegend(position = :rt)
	
	Permfig
end

# ╔═╡ fee7f74a-208e-48b7-adc1-5573b3d65dca
# ╠═╡ disabled = true
#=╠═╡
function copy_texmaker_format(t_range_optim_coupled,K_avg_norm_optim_coupled)
    
    lines = [@sprintf("%.6f %.6f \\\\", t_range_optim_coupled_pv[i], K_avg_norm_optim_coupled[i]) for i in 1:length(t_range_optim_coupled)]
    output_text = join(lines, "\n")
    
    clipboard(output_text)
    println("Copied to clipboard! Ready to paste into Texmaker.")
    println("First few lines:")
    println(output_text[1:min(100, length(output_text))])
end
  ╠═╡ =#

# ╔═╡ d9c906ea-022e-4d5e-a2cf-3490ca9950e6
# ╠═╡ disabled = true
#=╠═╡
copy_texmaker_format(t_range_optim_coupled,K_avg_norm_optim_coupled)
  ╠═╡ =#

# ╔═╡ Cell order:
# ╠═3c915980-5601-11f0-34f2-37291524e5a4
# ╠═ca47908b-cfef-4ced-ab9f-3e8d89e69970
# ╠═8c4ee984-9339-4f9e-af82-66ae9c342ce6
# ╠═a417ee46-af67-4d9e-af56-d3bdb621a2c9
# ╠═53beaddf-502a-4464-8891-0a4928739e2f
# ╠═312fa5d6-4a9d-4154-95a8-66820e555fea
# ╠═750a7e52-9e4c-4eb9-9a68-622dca3ae870
# ╠═d883de23-c596-423b-8866-e457a232a8d9
# ╠═b2525aea-5c06-416a-9efe-dfb34bfbdda9
# ╠═f9dadcaa-72e8-4383-afc5-cca26b6a5673
# ╠═f5ab7022-a472-414c-99c0-9e9f0181382c
# ╠═12c3b492-8f63-4b3e-a1d7-d0f346865eea
# ╠═b11fc159-2eaf-46d4-ac65-c5a40ab2ed6a
# ╠═719d7bbf-d6b2-49e1-b5cf-f9ed59f40953
# ╠═a8130b51-b286-4768-96ff-949dcccc55df
# ╠═a87a36eb-fd73-464a-9177-3b6a18fd7759
# ╠═8a3cd7ea-29fc-46be-b74f-26a282166017
# ╠═3f436653-0912-423b-b087-7888ccbef80c
# ╠═357ea7cb-9e0a-4e6e-98f2-f124fcedf5fb
# ╠═469194e5-0e01-439a-9cca-8dcbec018c1f
# ╠═2452d35f-f2b7-4202-a989-decbc58b260d
# ╠═07bd9879-d8d8-4c27-94c4-a6ce4437ba55
# ╠═d4a9accb-d00d-459d-ac7b-e15ea50dd71d
# ╠═fee7f74a-208e-48b7-adc1-5573b3d65dca
# ╠═d9c906ea-022e-4d5e-a2cf-3490ca9950e6
