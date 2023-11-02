mutable struct VerletList{Dim}
    box_list::Vector{Int64}
    # box_list[i] is either another index which is in the same box as 'i' or < 1

    heads::Array{Int64,2}
    # containts the first index out of each box

    # domain
    bottomleft::SVector{Dim,Float64}
    topright::SVector{Dim,Float64}
    dx_inv::SVector{Dim,Float64}
    strides::SVector{Dim,Int64}
end


function VerletList{Dim}(bottomleft, topright, R) where Dim

    grid_size = floor.(Int64, (topright - bottomleft) ./ R)
    dx = (topright - bottomleft) ./ grid_size

    box_list = Array{Int64}(undef, 0)
    heads = zeros(Int64, Int64.(grid_size)...)  # Array{1,Int64}(undef, n_cells)

    strides = SVector{2,Int64}(cumprod(vcat(1, size(heads)[1:end-1]...)))

    VerletList{Dim}(box_list, heads,
                    SVector{Dim,Float64}(bottomleft), SVector{Dim,Float64}(topright),
                    1. ./ dx,
                    strides)
end


struct VerletPairs{Dim}
    verletlist::VerletList{Dim}
end

Base.@kwdef struct VerletPairsIterator{Dim}
    i::Int64
    j::Int64
    box_i::CartesianIndex{Dim}
    box_j::CartesianIndex{Dim}
    neighbour_ind::Int64
end

function Base.iterate(vp::VerletPairs{Dim}) where {Dim}
    box_i = findfirst(x -> x > 0, vp.verletlist.heads)
    if isnothing(box_i)
        return nothing
    else
        i = vp.verletlist.heads[box_i]
        vpi = VerletPairsIterator(i=i, j=i,
            box_i=box_i, box_j=box_i,
            neighbour_ind=1)
        return iterate(vp, vpi)
    end
end

const neighbours_2d = CartesianIndex.([(0,0),(-1,1),(0,1),(1,1),(1,0)])

function Base.iterate(vp::VerletPairs{Dim}, s::VerletPairsIterator{Dim}) where {Dim}#::Union{Nothing,Tuple{Tuple{Int64,Int64},VerletPairIterator}}

    # get next element
    j = vp.verletlist.box_list[s.j]
    i = s.i
    box_i = s.box_i
    box_j = s.box_j
    neighbour_ind = s.neighbour_ind

    # if j < 0, then the according box is empty and we have to go the the next box
    while j <= 0
        # increment to next non-empty neighbour box
        while j <= 0 && neighbour_ind < length(neighbours_2d)
            neighbour_ind += 1
            box_j = box_i + neighbours_2d[neighbour_ind]

            # check if the box is within the grid
            if  box_j[1] >= 1 &&
                box_j[1] <= size(vp.verletlist.heads, 1) &&
                box_j[2] <= size(vp.verletlist.heads, 2)

                # go to first element of that new box
                j = vp.verletlist.heads[box_j[1], box_j[2]]
            end
        end

        # if last neighbour box is reached, increase 'i'
        if j <= 0 && neighbour_ind >= length(neighbours_2d)

            # all neighbours are done, increment 'i'
            i = vp.verletlist.box_list[i]
            neighbour_ind = 1

            # this loops move through all CartesianIndices(vl.heads)
            # I use the manual loop to avoid allocations
            while i <= 0
                box_i = box_i + CartesianIndex{2}(1,0)
                if box_i[1] > size(vp.verletlist.heads,1)
                    box_i = CartesianIndex{2}(1,box_i[2]+1)
                    if box_i[2] > size(vp.verletlist.heads,2)
                        return nothing
                    end
                end
                i = vp.verletlist.heads[box_i[1],box_i[2]]
            end
            j = i
        end
    end

    if i == j
        # skip self interaction
        return iterate(vp, VerletPairsIterator(i, j,
            box_i, box_j,
            neighbour_ind))
    end

    return ((i,j),VerletPairsIterator(i, j,
        box_i, box_j,
        neighbour_ind))  # at this line allocations occour!
end

Base.IteratorSize(vp::VerletPairs) = Base.SizeUnknown()


function update_lists!(vl::VerletList{Dim}, X) where Dim
    N = size(X,2)


    if N != size(vl.box_list,2)
        resize!(vl.box_list, N)
    end
    vl.heads .= 0
    try
        for i = 1:N
            box_idx = 1
            for k = 1:Dim
                box_idx += floor(Int64,(X[k,i] - vl.bottomleft[k])*vl.dx_inv[k]) * vl.strides[k]
            end
            vl.box_list[i] = vl.heads[box_idx]
            vl.heads[box_idx] = i
        end
    catch e
        idx = findfirst( x -> any(x .< vl.bottomleft) || any(x .> vl.topright), [X[:,i] for i=1:N])
        if !isnothing(idx)
            error("VerletList: Point #$(idx) is not in verlet list domain. This might indicate a numerical instability.
            X = $(X[:,idx]); VerletList Domain: (bottomleft=$(vl.bottomleft), topright=$(vl.topright)).")
        else
            throw(e)
        end
    end
    nothing
end

function adapt_bounds!(vl::VerletList{Dim}, s) where {Dim}
    N = size(s.X,2)
    out_of_bounds = []

    center = 0.5*(vl.topright + vl.bottomleft)
    rect_size = (vl.topright - vl.bottomleft)
    for i = 1:N
        if any(s.X[:,i] .< center - 0.4*rect_size) || any(s.X[:,i] .> center + 0.4*rect_size)
            push!(out_of_bounds, i)
        end
    end

    if length(out_of_bounds) > 0
        # @warn "Needed to adapt VerletList Domain."

        rect = Rectangle{Dim}(rect_size, center)

        extend!(rect, s.X[:,out_of_bounds])
        rect.size = 1.2 .* rect.size
        vl.topright = rect.center + 0.5*rect.size
        vl.bottomleft = rect.center - 0.5*rect.size

        R_max = 2 * maximum( cell -> cell.R_soft, s.cells )
        grid_size = floor.(Int64, (vl.topright - vl.bottomleft) ./ R_max)
        dx = (vl.topright - vl.bottomleft) ./ grid_size
        vl.dx_inv = 1. ./ dx

        vl.box_list = Array{Int64}(undef, 0)
        vl.heads = zeros(Int64, Int64.(grid_size)...)  # Array{1,Int64}(undef, n_cells)

        vl.strides = SVector{2,Int64}(cumprod(vcat(1, size(vl.heads)[1:end-1]...)))
    end
end

function get_box(vl::VerletList, ind)
    box = [vl.heads[ind]]
    while box[end] != 0
        push!(box, vl.box_list[box[end]])
    end
    return box[1:end-1]
end


function get_boxes(vl::VerletList)
    return reshape([get_box(vl,ind) for ind in eachindex(vl.heads)], size(vl.heads))
end
