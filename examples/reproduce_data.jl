using SEMTor, DataFrames, ProgressMeter, CSV

n_rep = 100
p_inm = 30
p_run = 50

printstyled("Load inputs ", bold=true, color=:blue)
printstyled("  Parameters: n_rep = $n_rep, p_inm = $p_inm, p_run = $p_run\n", color=:blue)

p_xlsx = SEMTor.load_parameters("input/ensemble_parameters.xlsx")
setvalue!(p_xlsx, "stat.num_rep", n_rep)
setvalue!(p_xlsx, "exp.run.sim_time_start", "(0,0, $(100 - p_run)% Inf)")
setvalue!(p_xlsx, "exp.INM.sim_time_start", "(0,0, $(100 - p_inm)% Inf)")

p = StaticParameters(p_xlsx)

event_states = Vector{Any}(undef, n_rep)
df = Vector{DataFrame}(undef, n_rep)

printstyled("Simulating ", bold=true, color=:blue)

prog = Progress(n_rep,1,"Simulate ")
Threads.@threads for rep in 1:n_rep
    obs = []
    states = simulate(p, rep; inplace=false, obs = obs, obs_callback = collect_times)
    df[rep] = analyse_state(states, p, rep, cell_type, obs)
    next!(prog)
end
GC.gc()

df[!,:scenario] = categorise.(df.A, df.B, df.P, df.R)
scenarios = unique(df.scenario)
sort_scenarios!(scenarios)

Δt = 6.0
Δy = 0.33

