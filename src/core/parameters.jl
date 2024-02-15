Base.@kwdef struct EpitheliumParameters

    N_init::Int64 = 30

    basal_curvature::Float64 = 0.0
    basal_curvature_elliSEMTor_factor::Float64 = 4.0
    init_curve_start::Float64 = 0.0
    curve_growth::Float64 = 0.1
    max_cytoskeleton_length::Float64 = Inf64
    running_zone::SVector{2,Float64} = [-2.0, 0.0]
    init_method::CellInitTypes = random
    init_distr::String = "control: 1"
    loop::Int = 0


    init_apical_junction_dist::Float64 = 1. / 3. # initialization of apical points
    init_basal_junction_dist::Float64 = 1. / 3.  # initialization of basal points
    init_zone::Rectangle{2} = Rectangle{2}(;size=[N_init * 1. / 3., 8.], center=[0., 4.])

    max_basal_junction_dist::Float64 = 1. / 3.   # maximum distance due to lateral adhesion, basal adhesion, rel to R_avg
    prob_out_div::Float64 = 0.8

    mu::Float64 = 0.02
    basal_damping_ratio::Float64 = 1.  # same damping as everything else
end

Base.@kwdef struct PBD
    alg_name::String = "PBD"
    dt::Float64 = 0.05
    n_substeps::Int64 = 40
end

Base.@kwdef struct DAHA
    alpha::Float64 = 0.01^2
    beta::Float64 = 0.25
    gamma::Float64 = 0.1^2
    rho::Float64 = 0.25
    eps::Float64 = 1e-6 
    c::Float64 = 1.0
    dt::Float64 = 1.0
end

Base.@kwdef struct StatisticalParameters
    num_rep::Int64 = 50
    random_seed::Int64 = 0
end

Base.@kwdef struct SimulationParameters
    name::String = "simulation"
    dt::Float64 = 0.1
    t_end::Float64 = 48.
end

Base.@kwdef struct StaticParameters
    epi::EpitheliumParameters = EpitheliumParameters()
    sim::SimulationParameters = SimulationParameters()
    alg::PBD = PBD()
    daha::DAHA = DAHA()
    stat::StatisticalParameters = StatisticalParameters()
    cell_types::Vector{CellTypeSampler} = [CellTypeSampler(EpiCellType(),RandomCellTypeParam[])]
end

cell_type_names(p::StaticParameters) = [cts.prototype.name for cts in p.cell_types]
cell_type_sampler_idx(p::StaticParameters, name::String) = findfirst(s -> s == name, cell_type_names(p))
cell_type_sampler(p::StaticParameters, name::String) = p.cell_types[cell_type_sampler_idx(p,name)]
cell_type_sampler(p::StaticParameters, i::Int64) = p.cell_types[i]
generate_cell_type(p::StaticParameters, name::String) = cell_type_sampler(p,name)()

function tryparsenan(T::Type{Float64}, x)
    v = tryparse(T, x)
    isnothing(v) && return NaN64
    return v
end

function tryparsenan(T::Type{Int64}, x)
    v = tryparse(T, x)
    isnothing(v) && return -1
    return v
end

# conversation of Dict to AllParameters
construct(::Type{T}, x::T; kw...) where {T} = x
construct(::Type{Float64}, s::String; nan_on_error=false, kw...) = nan_on_error ? tryparsenan(Float64, s) : parse(Float64, s)
construct(::Type{Int64}, s::String; nan_on_error=false, kw...) = nan_on_error ? tryparsenan(Int64, s) : parse(Int64, s)
construct(::Type{Float64}, s::Int64; kw...) = Float64(s)
construct(::Type{Symbol}, s::String; kw...) = Symbol(s)
function construct(::Type{SVector{2,Float64}}, s::String; kw...)
    return SVector{2,Float64}( Meta.eval(Meta.parse(s)))
end

function construct(::Type{T}, s::String; kw...) where {T <: Union{CellInitTypes,ReferenceTimePoint}}
    inst = instances(T)
    enum_dict = Dict(zip(string.(inst), inst))
    if haskey(enum_dict, s)
        return enum_dict[s]
    else
        @error "Unrecogniced option '$(s)' for type '$(T)'. Possible values are $(string.(inst))."
    end
end

function construct(::Type{Function}, s::String; kw...)
    if s == "divide_cell"
        return divide_cell
    elseif s == "reset_cell_cycle"
        return reset_cell_cycle
    elseif s == "loss_basal_connection"
        return loss_basal_connection
    elseif s == "loss_apical_connection"
        return loss_apical_connection
    else
        @error "Cannot convert '$(s)' into a julia function.\n  Options are: [divide_cell, reset_cell_cycle, loss_apical_connection, loss_basal_connection]"
    end
end

function construct(::Type{StaticArrays.SVector{2, Float64}}, d::Dict; kw...)
    if haskey(d, "x") && haskey(d, "y")
        return @SVector( Float64[d["x"], d["y"]] )
    elseif haskey(d, "w") && haskey(d, "h")
        return @SVector( Float64[d["w"], d["h"]] )
    else
        @error "Cannot construct 2D vector from $(d)."
    end
end

function construct(::Type{SEMTor.AbstractAlgorithm}, x::Dict; kw...)
    if !haskey(x, "alg_name")
        @error "Cannot find value for 'alg.alg_name' to determine numerical method."
    end

    alg_name = x["alg_name"]
    if alg_name == "PBD"
        return construct(PBD, x; kw...)
    elseif alg_name == "PGM"
        return construct(PGM, x; kw...)
    else
        @error "Algorithm $(alg_name) is currently not supported."
    end
end

function construct(::Type{Vector{T}}, V::Vector{A}; kw...) where {T<:Union{CellEvent,SpecialCellEvent,CellTypeSampler}, A}
    return [construct(T, x; idx=i, kw...) for (i,x) in enumerate(V)]
end

dict_types = Union{
            StaticParameters, SimulationParameters, EpitheliumParameters, StatisticalParameters,
            TimeInterval, Rectangle{2}, CellEvent, SpecialCellEvent, EpiCellType, PBD, DAHA
            }

function construct(::Type{T}, x::Dict{K,V}; idx = 1, kw...) where {T <: dict_types, K<:Union{Any,Symbol,String}, V}
    res = Dict{Symbol,Any}()
    types = Dict( f => t for (f,t) in zip(fieldnames(T),fieldtypes(T)) )

    # get default parameters
    if length(methods(T,(Dict,))) > 0
        def = T()
        for k in fieldnames(T)
            res[k] = getfield(def, k)
        end
    end
    # find all other parameters from the dict
    for (k, v) in x
        ks = Symbol(k)
        if haskey(types, ks)
            VT = types[ks]
            res[ks] = construct(VT, v; kw...)
        elseif ks in (:cell_types,)
            # ignore this field
        else
            @warn "Ignore unknown property $(ks) (with value $(v)) for type $(T)."
        end
    end

    if T == EpiCellType
        res[:prototype_idx] = idx
    end
    T(;res...)
end

function construct(::Type{CellTypeSampler}, x::T; idx=1, kw...) where {T <: Dict}
    prototype = construct(EpiCellType, x; nan_on_error=true, idx, kw...)
    rps = get_random_params(prototype, x)

    return CellTypeSampler(prototype, rps)
end


StaticParameters(d::Dict{K,Any}) where {K <: Union{Any,String}} = construct(StaticParameters, d)

EpiCellType(d::Dict{String,Any}) = construct(EpiCellType, d)
CellEvent(d::Dict{String,Any}) = construct(CellEvent, d)
SpecialCellEvent(d::Dict{String,Any}) = construct(SpecialCellEvent, d)
AffineTraffo(d::Dict{String,Any}) = construct(AffineTraffo, d)
TimeInterval(d::Dict{String,Any}) = construct(TimeInterval, d)
