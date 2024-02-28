using SEMTor, DataFrames, ProgressMeter, CSV, Statistics, XLSX
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

n_rep = get(ARGS, 1, 15000)
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

df[!,:scenario] = categorise.(df.A, df.B, df.S)


Δt = 6.0
Δy = 0.33

CSV.write("ensemble_data.csv", df)
save_binned_data("extrusion_statistics", df; bin_pos_emt = true, Δt = Δt, Δy = Δy)


using CSV, DataFrames, Statistics

function cor_(X,Y)
    inds = @. isfinite(X) && isfinite(Y)
    X_ = X[ inds ]
    Y_ = Y[ inds ]
    
    return sum(inds) > 0 ? cor(X_, Y_) : NaN
end


function delta_emt(df) 
    minmax = filter(isfinite, (df.A, df.B, df.S))
    if length(minmax) <= 1
        return NaN
    else
        return maximum(minmax) - minimum(minmax)
    end
end


df_cor = deepcopy(df)

df_cor[!,:Δt_emt] = delta_emt.(eachrow(df))

# rename columns and make events B and P mutually exclusive

df_cor.A .= @. isfinite(df.A)
df_cor.B .= @. isfinite(df.B) && !isfinite(df.P)
df_cor.P .= @. isfinite(df.B) && isfinite(df.P)
df_cor.S .= @. isfinite(df.S)

df_cor[!,:t_A] .= df.A 
df_cor[!,:t_B] .= df.B
df_cor[!,:t_P] .= df.B 
df_cor[!,:t_S] .= df.S

# make events B and P mutually exclusive
df_cor.t_B[.!df_cor.B] .= NaN
df_cor.t_P[.!df_cor.P] .= NaN

df_cor.y_B[.!df_cor.B] .= NaN
df_cor.y_P .= df.y_B
df_cor.y_P[.!df_cor.P] .= NaN


df_inm = df_cor[df.INM .== true, :]
df_no_inm = df_cor[df.INM .== false, :]

using CairoMakie

xaxis = ("A", "B", "P", "S", "t_A", "t_B", "t_P", "t_S", "Δt_emt", "y_init", "y_emt", "y_A", "y_B", "y_P", "y_S")

ydata_inm_a = [cor_(df_inm[!,col], df_inm.apical_extr) for col in xaxis]
ydata_no_inm_a = [cor_(df_no_inm[!,col], df_no_inm.apical_extr) for col in xaxis]

ydata_inm = [cor_(df_inm[!,col], df_inm.basal_extr) for col in xaxis]
ydata_no_inm = [cor_(df_no_inm[!,col], df_no_inm.basal_extr) for col in xaxis]

cor(df_inm.t_B, df_inm.basal_extr)

begin
    fig = Figure(size = (1024, 512), fontsize = 18)
    ax = Axis(fig[1, 1], ylabel = "correlation factor", title = "Apical extrusion", 
    xticks = (1:15,collect(xaxis)), xticklabelrotation = pi/4)
    ylims!(ax, -1, 1)

    xdata = Float64.(1:15) 
    hlines!(ax, [0], color = :gray, linestyle = :dash, linewidth = 2)
    scatter!(ax, ydata_inm_a, color = "#d6d5d3", markersize = 15, strokewidth = 1, strokecolor = :black, label = "INM")
    scatter!(ax, ydata_no_inm_a, color = "#ffa6a4", markersize = 15, strokewidth = 1, strokecolor = :black, label = "No INM")
    

    ax = Axis(fig[1, 2], title = "Basal extrusion", 
    xticks = (1:15,collect(xaxis)), xticklabelrotation = pi/4)
    ylims!(ax, -1, 1)
    hlines!(ax, [0], color = :gray, linestyle = :dash, linewidth = 2)
    scatter!(ax, ydata_inm, color = "#d6d5d3", markersize = 15, strokewidth = 1, strokecolor = :black, label = "INM")
    scatter!(ax, ydata_no_inm, color = "#ffa6a4", markersize = 15, strokewidth = 1, strokecolor = :black, label = "no INM")
    
    Legend(fig[1,end+1], ax, tellheight = false)

    save("correlation_factors.png", fig)
    fig
end