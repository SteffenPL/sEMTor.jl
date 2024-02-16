EpiCellVector = typeof(StructVector(EpiCell[]))

Base.@kwdef mutable struct State{Dim}
    X::HArray{Float64,Dim} = hzeros(Val(Dim),0)  # nuclei positions
    A::HArray{Float64,Dim} = hzeros(Val(Dim),0)  # apical positions
    B::HArray{Float64,Dim} = hzeros(Val(Dim),0)  # basal positions

    apicalcons::EdgeList{Float64} = EdgeList{Float64}()
    basalcons::EdgeList{Missing} = EdgeList{Missing}() # no custom data needed

    cells::EpiCellVector = StructVector(EpiCell[])  # array of epithelial cells

    t::Float64 = 0.  # time
    types::Vector{EpiCellType} = []  # array with the concrete instances of each cell type
end

const State2D = State{2}

function by_label(e, label)
    idx = findfirst(l -> l == label, e.cells.label)
    if isnothing(idx)
        return nothing
    else
        return e[idx]
    end
end

Dim(::State{n}) where {n} = n
is3D(s::State2D) = false

num_cells(s::State) = length(s.cells.R_soft)
get_type(s,i::Int64) = s.types[s[i].type_idx]
get_type(s,c::EpiCell) = s.types[c.type_idx]
get_type(s,c::LazyRow{EpiCell}) = s.types[c.type_idx]
get_prototype(s,p,i) = p.cell_types[s[i].prototype_idx]
get_age(s,i::Int64) = s.cells.age[i]
new_idx(s::State) = num_cells(s) + 1
new_label(s::State) = maximum(s.cells.label) + 1

function quickcopy(s::State{Dim}) where {Dim}
    return State{Dim}(
        X = deepcopy(s.X),
        A = deepcopy(s.A),
        B = deepcopy(s.B),
        apicalcons = deepcopy(s.apicalcons),
        basalcons = deepcopy(s.basalcons),
        cells = deepcopy(s.cells),
        t = s.t,
        types = s.types
    )
end

function new_cell(s::State,p::StaticParameters,i::Int64)
    p_idx = s.cells.prototype_idx[i]
    type = cell_type_sampler(p,p_idx)()
    push!(s.types, type)
    type_idx = length(s.types)
    EpiCell(type; label=new_label(s), type_idx)
end

function new_cell_copy_startimes(s::State,p::StaticParameters,i::Int64)
    p_idx = s.cells.prototype_idx[i]
    type = p.cell_types[p_idx](get_type(s,i))
    push!(s.types, type)
    type_idx = length(s.types)
    EpiCell(type;
        zip(fieldnames(EpiCell),
            getproperty(s[i],prop) for prop in fieldnames(EpiCell))...,
        label=new_label(s), type_idx)
end

function is_below_basal_layer(s::State{Dim}, p::StaticParameters, i::Int64, basal_info) where {Dim}
    epi = p.epi
    (;no_curvature, Rb, y0, Rc) = basal_info

    if no_curvature
        return s.B[2,i] < 0.0 - 0.01
    else
        if Rb > 0
            return sum(x->x^2, s.B[:,i] - Rc) > Rb^2 + 0.01
        else
            return sum(x->x^2, s.B[:,i] - Rc) < Rb^2 - 0.01
        end
    end
end

function nuclei_is_below_basal_layer(s::State{Dim}, p::StaticParameters, i::Int64, basal_info) where {Dim}
    epi = p.epi
    (;no_curvature, Rb, y0, Rc) = basal_info

    if no_curvature
        return s.X[2,i] < 0.0 - 0.01
    else
        if Rb > 0
            return sum(x->x^2, s.X[:,i] - Rc) > Rb^2 + 0.01
        else
            return sum(x->x^2, s.X[:,i] - Rc) < Rb^2 - 0.01
        end
    end
end

function get_cell_event(s, i, event_name)
    t = get_type(s, i)
    for e in t.cell_events
        if e.name == event_name
            return e
        end
    end
    for e in t.special_cell_events
        if e.name == event_name
            return e
        end
    end
    return nothing
end

import Base.getindex
import Base.setindex!
import Base.firstindex
import Base.lastindex
import Base.length
import Base.iterate


getindex(state::State, idx::Int64) = LazyRow(state.cells, idx)
getindex(state::State, I::UnitRange{Int64}) = map(i -> state[i], I)
getindex(state::State, ::Colon) = state.cells

getindex(state::State, I::UnitRange{Int64}, sym::Symbol) = map( i -> getproperty(state[i], sym), I)
getindex(state::State, ::Colon, sym::Symbol) = getproperty.(state.cells, sym)

function setindex!(state::State, cell::EpiCell, idx::Int64)
    for prop in fieldnames(EpiCell)
        getproperty(state.cells, prop)[idx] = getproperty(cell, prop)
    end
end

firstindex(s::State) = 1
lastindex(s::State) = num_cells(s)
length(s::State) = num_cells(s)

function iterate(s::State)
    if length(s) > 1
        return s[1], 1
    else
        return nothing
    end
end

function iterate(s::State, iter::Int64)
    if length(s) > iter
        return s[iter+1], iter+1
    else
        return nothing
    end
end

function Base.show(io::IO, state::State)
    compact = get(io, :compact, false)
    if !compact
        @printf(io, "State%dD at t=%.2f [h] with %d cells.", Dim(state), state.t[], length(state))
    else
        @printf(io, "State%dD(t=%.1f,%d)", Dim(state), state.t[], length(state))
    end
end


function (states::Vector{<:State})(t::Float64)
    ind = findfirst( s -> s.t >= t - 0.01, states )
    return isnothing(ind) ? last(states) : states[ind]
end