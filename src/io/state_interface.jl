function getproperty_save(x, sym::Symbol)
    if hasproperty(x, sym)
        return getproperty(x, sym)
    else
        return missing
    end
end

function create_cell_dataframe(state::State2D)
    cs = state.cells
    DataFrame(  label = state[:,:label],
                idx = collect(eachindex(cs)),
                X1 = state.X[1,:], X2 = state.X[2,:],
                A1 = state.A[1,:], A2 = state.A[2,:],
                B1 = state.B[1,:], B2 = state.B[2,:],
                R_soft = cs.R_soft,
                R_hard = cs.R_hard,
                type =      map( cell -> get_type(state, cell).name,                state.cells),
                age = cs.age,
                time_division = state.cells.division_time,
                has_apical_cytoskeleton = cs.has_apical_cytos,
                apical_rest_length = cs.apical_cytos_rest_length,
                apical_cytos_strain = cs.apical_cytos_strain,
                has_basal_cytoskeleton = cs.has_basal_cytos,
                basal_rest_length = cs.basal_cytos_rest_length,
                basal_cytos_strain = cs.basal_cytos_strain,
                k_cytoskeleton = cs.k_cytoskeleton,
                k_apical_junction = cs.k_apical_junction,
                stiffness_repulsion = cs.stiffness_repulsion,
                stiffness_nuclei_apical = cs.stiffness_nuclei_apical,
                stiffness_nuclei_basal = cs.stiffness_nuclei_basal,
                stiffness_apical_apical = cs.stiffness_apical_apical,
                stiffness_straightness = cs.stiffness_straightness,
                constrained_to_basal_layer = cs.on_basal_layer,
                diffusion = cs.diffusion,
                duration_G2 = cs.duration_G2,
                duration_mitosis = cs.duration_mitosis,
                time = fill(state.t, length(state))
                )
end

function create_connection_dataframe(state::State2D)
    cons = collect(bio_connections(state))
    df = DataFrame( type = map( con -> string(typeof(con)), cons),
                u = map( con -> state[con.u.idx].label, cons),
                v = map( con -> state[con.v.idx].label, cons),
                rest_length =           map( con -> getproperty_save(con, :rest_length),         cons),
                desired_rest_length =   map( con -> getproperty_save(con, :desired_rest_length), cons),
                time = fill(state.t, length(cons))
             )
    categorical!(df, :type)
    sort!(df,[:type])
    df.type = string.(df.type)
    return df
end


function create_cell_dataframe(states::Array{State2D,1}; rep::Union{Int64,Missing}=missing)
    df_iter = DataFrame()
    if !ismissing(rep)
        df_iter[!, :rep] = fill(rep,length(states))
    end
    df_iter[!,:iter] = 1:length(states)
    df_iter[!,:time] = map(state -> state.t, states)

    join(df_iter, vcat(map(state -> create_cell_dataframe(state), states)...), on=:time)
end

function create_cell_dataframe(sim_reps::Array{Array{State2D,1},1})
    vcat((create_cell_dataframe(state, rep=rep) for (rep,state) in enumerate(sim_reps))...)
end

function create_connection_dataframe(states::Array{State2D,1}; rep::Union{Int64,Missing}=missing)
    df_iter = DataFrame()
    if !ismissing(rep)
        df_iter[!, :rep] = fill(rep,length(states))
    end
    df_iter[!,:iter] = 1:length(states)
    df_iter[!,:time] = map(state -> state.t, states)

    join(df_iter, vcat(map(state -> create_connection_dataframe(state), states)...), on=:time)
end

function create_connection_dataframe(sim_reps::Array{Array{State2D,1},1})
    vcat((create_connection_dataframe(state, rep=rep) for (rep,state) in enumerate(sim_reps))...)
end
