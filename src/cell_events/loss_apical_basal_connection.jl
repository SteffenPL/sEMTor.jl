function loss_apical_connection(state::State, p::StaticParameters, i::Int64)
    cell = state[i]

    # we could also deactive the cytoskeleton alltogether, but we
    # instead let it degenerate:
    #cell.has_apical_cytos = false
    cell.apical_cytos_strain = -1.
    cell.stiffness_nuclei_apical *= 0.1

    adjs = findall( edge -> i in edge, state.apicalcons.edges)

    m = length(adjs)
    @assert m <= 2 """In 2D, we only expect maximal 2 apical neighbours. Cell #$(i), length(adjs) = $m.
    t = $(state.t), i = $i, adjs = $adjs"""

    edges = state.apicalcons.edges[adjs]
    rest_lengths = state.apicalcons.data[adjs]
    if m == 2
        # remove one edge and add the rest lengths of the other edge
        if edges[1][1] == i
            # v => i => u
            v = edges[2][1]
            u = edges[1][2]
            state.apicalcons.edges[adjs[1]] = v => u
        else # edges[1][2] == i
            # v => i => u
            v = edges[1][1]
            u = edges[2][2]
            state.apicalcons.edges[adjs[1]] = v => u
        end
        rest_length = rest_lengths[1] + rest_lengths[2]
        state.apicalcons.data[adjs[1]] = rest_length


        deleteat!( state.apicalcons.edges, adjs[2])
        deleteat!( state.apicalcons.data, adjs[2])
    elseif m == 1
        deleteat!( state.apicalcons.edges, adjs[1])
        deleteat!( state.apicalcons.data, adjs[1])
    end

end


function loss_basal_connection(state::State, p::StaticParameters, i::Int64)
    cell = state[i]

    # we could also deactive the cytoskeleton alltogether, but we
    # instead let it degenerate:
    #cell.has_apical_cytos = false
    cell.basal_cytos_strain = -1.
    cell.stiffness_nuclei_basal *= 0.1

    adjs = findall( edge -> i in edge, state.basalcons.edges)

    m = length(adjs)
    @assert m <= 2 "In 2D, we only expect maximal 2 apical neighbours."

    edges = state.basalcons.edges[adjs]
    if m == 2
        if edges[1][1] == i
            # v => i => u
            u = edges[1][2]
            v = edges[2][1]
            state.basalcons.edges[adjs[1]] = v => u
        else # edges[1][2] == i
            # v => i => u
            v = edges[1][1]
            u = edges[2][2]
            state.basalcons.edges[adjs[1]] = v => u
        end
        deleteat!( state.basalcons.edges, adjs[2])
        deleteat!( state.basalcons.data, adjs[2])
    elseif m == 1
        deleteat!( state.basalcons.edges, adjs[1])
        deleteat!( state.basalcons.data, adjs[1])
    end
end
