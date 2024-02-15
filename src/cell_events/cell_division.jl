function divide_cell(state::State{d}, p::StaticParameters, i::Int64) where {d}

    cell = state[i]

    # store value which we want to keep
    acrl = cell.apical_cytos_rest_length
    bcrl = cell.basal_cytos_rest_length
    hac = cell.has_apical_cytos
    hbc = cell.has_basal_cytos

    cell1 = new_cell_copy_startimes(state,p,i)
    init_cell!(cell1, get_type(state,cell1), t=state.t)
    new_type = get_type(state, cell1)
    state[i] = cell1

    children = [cell1]
    childrenindices = [i]

    if rand() >= p.epi.prob_out_div
        cell2 = new_cell(state,p,i)
        push!(state.cells, cell2)
        push!(children, cell2)

        init_cell!(cell2, get_type(state,cell2), t=state.t)

        x = state.X[:,i]
        a = state.A[:,i]
        b = state.B[:,i]

        if d == 2
            div_dir = [1., 0.]
            offset = 0.05 * cell.R_soft * div_dir
        else
            div_dir_xz = sincos( rand() * 2*pi)
            div_dir = [div_dir_xz[1], 0., div_dir_xz[2]]
            offset = 0.05 * cell.R_soft .* div_dir
        end
    end

    for child in children
        child.apical_cytos_rest_length = acrl
        child.basal_cytos_rest_length = bcrl
        child.has_apical_cytos = hac
        child.has_basal_cytos = hbc
    end

    if length(children) > 1
        j = num_cells(state)
        push!(childrenindices, j)

        for k in 1:length(state.apicalcons)
            (u,v) = state.apicalcons.edges[k]
            if i == v
                state.apicalcons.edges[k] = u => i
            end
            if i == u
                state.apicalcons.edges[k] = j => v
            end
        end

        for k in 1:length(state.basalcons)
            (u,v) = state.basalcons.edges[k]
            if i == v
                state.basalcons.edges[k] = u => i
            end

            if i == u
                state.basalcons.edges[k] = j => v
            end
        end
        add_edge!(state.apicalcons, i, j, get_type(state,i).apical_junction_init)
        add_edge!(state.basalcons, i, j)

        state.X = cat(state.X, x+offset, dims=2)
        state.A = cat(state.A, a+offset, dims=2)
        state.B = cat(state.B, b+offset, dims=2)

        state.X[:,i] -= offset
        state.A[:,i] -= offset
        state.B[:,i] -= offset
    end


    for (i, child) in zip(childrenindices, children)
        for event in get_type(state,i).cell_events
            init_cell_event!(state, p, i, event)
        end
        for event in get_type(state,i).special_cell_events
            init_cell_event!(state, p, i, event)
        end
        state[i] = child
    end



    return false
end
