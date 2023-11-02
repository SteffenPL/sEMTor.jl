abstract type AbstractAlgorithm end
abstract type AbstractCellEvent end

@enum CellPhase G0=0 G1=1 G2=2 mitosis=3
@enum CellInitTypes random=0 clustered=1

HArray{T,Dim} = HybridArray{Tuple{Dim,StaticArrays.Dynamic()},T,2,2,Array{T,2}} where {T,Dim}
hzeros(::Val{Dim}, n) where {Dim} = HArray{Float64,Dim}(zeros(Dim,n))

Base.@kwdef mutable struct TimeInterval
    min::Float64 = -Inf64
    max::Float64 = Inf64
end

function guess_value(s)
    x = tryparse(Int64, s)
    !isnothing(x) && return x
    x = tryparse(Float64, s)
    !isnothing(x) && return x
    return nothing
end

Base.@kwdef mutable struct Rectangle{Dim}
    size::SVector{Dim,Float64}
    center::SVector{Dim,Float64}
end

function Rectangle{2}(;size = [1., 1.], center = [0.,0.])
    return Rectangle{2}(SVector{2,Float64}(size), SVector{2,Float64}(center))
end

function lims_to_rect(xlims = [-0.5, 0.5], ylims = [-0.5,0.5])
    return Rectangle{2}(@SVector[xlims[2]-xlims[1], ylims[2] - ylims[1]], @SVector[0.5*xlims[1]+0.5*xlims[2], 0.5*ylims[1] + 0.5*ylims[2]])
end

function Rectangle{3}(;size = [1., 1., 1.], center = [0.,0., 0.])
    return Rectangle{3}(SVector{2,Float64}(size), SVector{2,Float64}(center))
end

dimension(rect::Rectangle{2}) = 2
dimension(rect::Rectangle{3}) = 3

function transform_to_unit_rect(r::Rectangle, p)
    return (p .- r.center + 0.5 .* r.size) ./ r.size
end

function map_to_rect(r::Rectangle, p)
    return r.size .* p - 0.5 .* r.center
end

function extend!(r::Rectangle, p)
    a = minimum( hcat(minimum(p, dims=2)[:], r.center - 0.5*r.size), dims=2)[:]
    b = maximum( hcat(maximum(p, dims=2)[:], r.center + 0.5*r.size), dims=2)[:]
    r.center = 0.5*(b+a)
    r.size = (b-a)
    nothing
end

function extend!(r::Rectangle, r2::Rectangle)
    extend!(r, [r2.center - 0.5*r2.size, r2.center + 0.5*r2.size])
end

get_lims(r::Rectangle, dim::Int64) = r.center[dim] .+ r.size[dim] .* (-0.5,0.5)
get_xlims(r::Rectangle) = get_lims(r,1)
get_ylims(r::Rectangle) = get_lims(r,2)
get_zlims(r::Rectangle) = get_lims(r,3)

sqdist(x,y) = sum(z -> z^2, x-y)
dist(x,y) = sqrt(sum(z -> z^2, x-y))
dist(x) = sqrt(sum(z -> z^2, x))
