struct EdgeList{TData}
    edges::Vector{Pair{Int64, Int64}}
    data::Vector{TData}
end

EdgeList{T}() where {T} = EdgeList{T}([],Array{T}(undef,0))

function Base.empty!(el::EdgeList)
    empty!(el.edges)
    empty!(el.data)
    nothing
end

Base.length(el::EdgeList) = length(el.edges)

@inline EdgeIter(g::EdgeList{TData}) where {TData} = zip(g.edges, g.data)

edge(u,v) = Pair(minmax(u,v)...)

function add_edge!(g::EdgeList{T}, u::Int64, v::Int64, data::Union{T,Missing}=missing) where {T}
    push!(g.edges, edge(u,v))
    push!(g.data, data)
end

function delete_edge!(g::EdgeList, u::Int64, v::Int64)
    deleteidx = edge(u,v) .== g.edges
    deleteat!(g.edges, deleteidx)
    deleteat!(g.data, deleteidx)
end

function edges_with_data(g::EdgeList{T}) where {T}
    return zip(g.edges, g.data)
end

function getdata(g::EdgeList{T}, u::Int64, v::Int64) where {T}
    a = edge(u,v)
    idx = findfirst(b -> a == b, g.edges)
    return g.data[idx]
end

function setdata(g::EdgeList{T}, u::Int64, v::Int64, d::T) where {T}
    a = edge(u,v)
    idx = findfirst(b -> a == b, g.edges)
    g.data[idx] = d
end
