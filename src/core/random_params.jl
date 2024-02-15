struct RandomCellTypeParam
    des_param::Tuple{Symbol,Int64,Symbol}
    des_param_str::String
    src_val::Float64
    src_param::Tuple{Symbol,Int64,Symbol}
    src_param_str::String
    src_rand::Tuple{Symbol,Float64,Float64,Float64,Float64}
end

function Base.show(io::IO, rct::RandomCellTypeParam)
    print(io, string(rct.src_param_str, "→", rct.des_param_str, " = ", rct.src_val))
end

function generate_random_param!(ct, rp::RandomCellTypeParam)
    if isnan(ct[rp.des_param])
        val = rp.src_val
        if rp.src_param[1] != :nothing && rp.src_rand[1] != :uniform_start_dep
            val += ct[rp.src_param]
        end
        if rp.src_rand[1] == :uniform
            val += rand(Uniform(rp.src_rand[2],rp.src_rand[3]))
        end
        if rp.src_rand[1] == :uniform_start_dep && rp.src_param[1] != :nothing
            val += rand(Uniform(ct[rp.src_param],rp.src_rand[3]))
        end
        if rp.src_rand[1] == :uniform_inf
            if rand() < rp.src_rand[4]
                val = rp.src_rand[5]
            else
                if rp.src_rand[2] != rp.src_rand[3]
                    val += rand(Uniform(rp.src_rand[2],rp.src_rand[3]))
                else
                    val += rp.src_rand[2]
                end
            end
        end
        ct[rp.des_param] = val
    end
end

function RandomCellTypeParam(ct, des_param_str::String, val::String)
    des_param = to_param_triple(ct, des_param_str)
    @assert typeof(des_param) == Tuple{Symbol,Int64,Symbol} """ Cannot parse random parameter '$des_param_str = $val' for cell type $(ct.name)."""

    src_val = 0.0
    src_rand = (:nothing,0.0,0.0,0.0,0.0)
    src_param = (:nothing,0,:nothing)
    src_param_str = ""

    parts = strip.(split(val,'+'))
    for part in parts
        v = guess_value(part)
        if !isnothing(v)
            src_val += v
        elseif part[1] == '(' && part[end] == ')'
            vals = strip.(split(part[2:end-1],','))
            if length(vals) == 2
                v1 = guess_value(vals[1])
                v2 = guess_value(vals[2])

                if isnothing(v1) 
                    src_param_str = vals[1]
                    src_param = to_param_triple(ct, vals[1])
                    src_rand = (:uniform_start_dep,0., v2, 0.0, 0.0)
                else
                    src_rand = (:uniform,minmax(v1,v2)...,0.0,0.0)
                end
            elseif length(vals) == 3
                if '%' in vals[3]
                    p_vals = strip.(split(vals[3],'%'))
                    v1 = guess_value(vals[1])
                    v2 = guess_value(vals[2])
                    p3 = guess_value(p_vals[1]) / 100.0
                    v3 = guess_value(p_vals[2])
                    src_rand = (:uniform_inf,minmax(v1,v2)...,p3,v3)
                else
                    @error "Cannot parse '$val' since percentage of third value is missing. e.g. (6,12,30%Inf)"
                end
            else
                @error "Cannot parse '$val' as a random parameter."
            end
        else
            src_param_str = part
            src_param = to_param_triple(ct, part)
        end
    end
    return RandomCellTypeParam(des_param, des_param_str, src_val, src_param, src_param_str, src_rand)
end

function generate_random_celltype(prototype, random_params::Vector{RandomCellTypeParam})
    ct = deepcopy(prototype)
    for rp in random_params
        generate_random_param!(ct, rp)
    end
    return ct
end



function inherit_random_values!(new_type, old_type)

    for (idx,event) in enumerate(new_type.cell_events)
        if (event.cell_ref_time == Cell_Birth
            && event.cell_time_start == 0.0)
            event.sim_time_start = old_type.cell_events[idx].sim_time_start
        end
    end

    for (idx,event) in enumerate(new_type.special_cell_events)
        if (event.cell_ref_time == Cell_Birth
            && event.cell_time_start == 0.0)
            event.sim_time_start = old_type.special_cell_events[idx].sim_time_start
        end
    end
end

function generate_random_celltype(prototype, random_params::Vector{RandomCellTypeParam}, parent::EpiCellType)
    ct = deepcopy(prototype)
    inherit_random_values!(ct, parent)
    for rp in random_params
        generate_random_param!(ct, rp)
    end
    return ct
end



struct CellTypeSampler
    prototype::EpiCellType
    random_params::Vector{RandomCellTypeParam}
end


function Base.show(io::IO, cts::CellTypeSampler)
    compact = get(io, :compact, false)
    if !compact
        println(io, string("Random ", cts.prototype.name, " type:"))
        print(io, cts.prototype)
    else
        print(io, string("Random ", cts.prototype.name, " type"))
    end
end

function (cts::CellTypeSampler)()
    return generate_random_celltype(cts.prototype,cts.random_params)
end


function (cts::CellTypeSampler)(parent)
    return generate_random_celltype(cts.prototype,cts.random_params, parent)
end


is_random(cts::CellTypeSampler) = !isempty(random_params)

# getindex/setindex! for EpiCellType are basically pointers on the cell parameter
# which we want to change. Not pretty, but it works
function Base.getindex(e::EpiCellType, i::Tuple{Symbol,Int64,Symbol})
    if i[1] == :cell_events
        return getproperty(e.cell_events[i[2]],i[3])
    elseif i[1] == :special_cell_events
        return getproperty(e.special_cell_events[i[2]],i[3])
    elseif i[1] == :life_span
        if i[3] == :min
            return e.life_span.min
        else
            return e.life_span.min
        end
    else
        return getproperty(e,i[1])
    end
end

function Base.setindex!(e::EpiCellType, val, i::Tuple{Symbol,Int64,Symbol})
    if i[1] == :cell_events
        return setproperty!(e.cell_events[i[2]],i[3], val)
    elseif i[1] == :special_cell_events
        return setproperty!(e.special_cell_events[i[2]],i[3], val)
    elseif i[1] == :life_span
        if i[3] == :min
            e.life_span.min = val
        else
            e.life_span.min = val
        end
    else
        setproperty!(e,i[1], val)
    end
end

function to_param_triple(ct::EpiCellType, s)
    x = split(s, '.')
    length(x) == 1 && return (Symbol(x[1]),0,:nothing)
    y = x[((x[1] == ct.name) ? 2 : 1):end]

    if y[1] == "life_span"
        y[2] == "min" && return (:life_span,0,:min)
        y[2] == "max" && return (:life_span,0,:max)
    end

    for sym_ce in (:cell_events, :special_cell_events)
        i_ce = findfirst(ce -> ce.name == y[1], getproperty(ct,sym_ce))
        if !isnothing(i_ce)
            s = Symbol(y[2])
            if !hasfield(eltype(getproperty(ct,sym_ce)),s)
                @error """
                    Cannot determine random or dependent parameters.
                    Unknown cell event parameter '$(y[2])'.
                    Options are $(string.(fieldnames(SEMTor.CellEvent)))."""
            end
            return (sym_ce,i_ce,s)
        end
    end
    @error """
        Cannot determine random or dependent parameters.
        Unknown cell type parameter '$(join(x, ", "))'."""
    return nothing
end

# Parsing part
function same_ending(endings, s::String)
    for ending in endings
        n = length(ending)
        if length(s) >= n && (@view s[end-n+1:end]) == ending
            return true
        end
    end
    return false
end


function get_unknown_parameters(d)
    unknowns = Dict{String,String}()
    knowns = String[]

    # find all undetermined parameters
    for (key, val) in d
        if (val isa String &&
            isnothing(guess_value(val)) &&
            !same_ending(
            ("name","symbol","julia_function","cell_ref_time",
            "color","init_distr","init_method"),
            key))
            unknowns[key] = val
        else
            push!(knowns, key)
        end
    end

    return unknowns, knowns
end


function get_random_params(prototype::EpiCellType, d)
    unknowns, knowns = get_unknown_parameters(flatten(d))
    rps = [RandomCellTypeParam(prototype, k, v) for (k,v) in unknowns]
    # sort according to their dependencies
    n = length(rps)

    rps_ordered = RandomCellTypeParam[]
    for _ = 1:n
        for (k,rp) in enumerate(rps)
            if rp.src_param_str in knowns || rp.src_param[1] == :nothing
                push!(knowns, rp.des_param_str)
                push!(rps_ordered, rp)
                deleteat!(rps, k)
                break
            end
        end
    end

    @assert length(rps) == 0 """Could not resolve all dependencies of parameters for cell type '$(prototype.name)'.
    Missing values are '$(join([rp.des_param_str for rp in rps],", "))'.
    """

    return rps_ordered
end
