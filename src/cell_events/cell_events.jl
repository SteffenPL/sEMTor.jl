
@enum ReferenceTimePoint begin
    Cell_Birth
    Cell_G2_Start
    Cell_Mitosis_Start
    Cell_Division
end

Base.@kwdef mutable struct CellEvent <: AbstractCellEvent
    name::String = "R_hard"
    cell_ref_time::ReferenceTimePoint = Cell_Birth
    cell_time_start::Float64 = 0.
    sim_time_start::Float64 = 0.
    sim_time_end::Float64 = Inf64
    symbol::Symbol = :R_hard
    abs_value::Float64 = 0.0
    factor::Float64 = 1.0
end

function CellEvent(name::String, param_sym::Symbol, abs_value::Float64, factor::Float64 = 0.; kwargs...)
    @assert param_sym in fieldnames(EpiCell) "symbol $(param_sym) is not a cell parameter\n Options are:\n $(fieldnames(EpiCell))."
    CellEvent(;name=name, symbol=param_sym, abs_value, factor, kwargs... )
end

function CellEvent(param_sym::Symbol, abs_value::Float64, factor::Float64 = 0.; kwargs...)
    CellEvent(string(param_sym), param_sym, abs_value, factor; kwargs... )
end

no_event(s, ::Int64) = nothing

Base.@kwdef mutable struct SpecialCellEvent <: AbstractCellEvent
    julia_function::Function = no_event
    cell_ref_time::ReferenceTimePoint = Cell_Birth
    cell_time_start::Float64 = 0.
    sim_time_start::Float64 = 0.
    sim_time_end::Float64 = Inf64
    name::String = replace(string(julia_function), "SEMTor." => "")
end

SpecialCellEvent(f::Function; kwargs...) = SpecialCellEvent(;julia_function=f, kwargs...)


function Base.show(io::IO, e::CellEvent)
    print(io, e.name, ": Start: ", e.sim_time_start,)
end

function Base.show(io::IO, e::SpecialCellEvent)
    print(io, e.name, ": Start: ", e.sim_time_start,)
end