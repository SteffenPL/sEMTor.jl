using SEMTor, DataFrames, ProgressMeter, CSV, Statistics
include("analyze_state.jl")

#=  
    Ensebmle simulation of the model for reproducing the data in 

        S. Plunder, C. Danesin, B. Glise, M. A. Ferreira, S. Merino-Aceituno, E. Theveneau, 
        Modelling variability and heterogeneity of EMT scenarios highlights nuclear positioning and protrusions as main drivers of extrusion. 
        (2023)

    Figure 7.

    The script runs an ensemble simulation of the model where the time of EMT specific events, 
    the occurance of INM and the protrusive activity are randomly varied for EMT cells. 

    The inputs are: 
        - n_rep: number of simulations
        - p_inm: percentage of INM occurance
        - p_run: percentage of protrusive activity


    Notice that this code can be executed in parallel by using the command 
        julia -p n reproduce_data.jl 1000 50 30
    where n is the number of threads to use.

    Moreover, it is shared as an executable code ocean capsule: URL
=#

n_rep = get(ARGS, 1, 150000)
p_inm = get(ARGS, 2, 30)
p_run = get(ARGS, 3, 50)


# load parameter files and change random probabilities:
printstyled("Load inputs ", bold=true, color=:blue)
printstyled("  Parameters: n_rep = $n_rep, p_inm = $p_inm, p_run = $p_run\n", color=:blue)
p_xlsx = SEMTor.load_parameters("ensemble_parameters.xlsx")
setvalue!(p_xlsx, "stat.num_rep", n_rep)
setvalue!(p_xlsx, "exp.run.sim_time_start", "(0,0, $(100 - p_run)% Inf)")
setvalue!(p_xlsx, "exp.INM.sim_time_start", "(0,0, $(100 - p_inm)% Inf)")

# convert to parameter object
p = StaticParameters(p_xlsx)

# store the result in a dataframes 
df = Vector{DataFrame}(undef, n_rep)
cell_type = "exp"  # which cell type to analyse

printstyled("Simulating ", bold=true, color=:blue)

prog = Progress(n_rep,1,"Simulate ")
Threads.@threads for rep in 1:n_rep
    obs = []
    states = simulate(p, rep; inplace=false, obs = obs, obs_callback = collect_times)
    df[rep] = analyse_state(states, p, rep, cell_type, obs)
    next!(prog)

    if rep % 1000 == 0
        GC.gc()  # just to be sure
    end
end

df = reduce(vcat, df)

df[!,:scenario] = categorise.(df.A, df.B, df.P, df.R)
scenarios = unique(df.scenario)
sort_scenarios!(scenarios)

Δt = 6.0
Δy = 0.33

CSV.write("ensemble_data.csv", df)

