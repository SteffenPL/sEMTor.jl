function init_state(p::StaticParameters, rep = 1)
    Random.seed!(p.stat.random_seed + rep - 1)
    X, A, B = init_positions(p)

    init_cell_type_names = generate_init_cell_types(p)
    all_cell_types = [cts.prototype.name for cts in p.cell_types]
    N = p.epi.N_init

    types = EpiCellType[]
    cells = StructVector(EpiCell[])

    for idx in 1:N
        if !(init_cell_type_names[idx] in all_cell_types)
            @error """Unknown type $(init_cell_type_names[idx]) in 'init_distr'. Could it be that 'init_distr' is not correct? \n
            ('init_distr' is $(p.epi.init_distr); known cell types are $(all_cell_types).)"""
        end
        new_type = generate_cell_type(p, init_cell_type_names[idx])
        push!(types, new_type)

        type_idx = length(types)
        epi_cell = EpiCell(new_type; label = idx, type_idx=type_idx)
        init_cell!(epi_cell, new_type)
        push!(cells, epi_cell)
    end

    s = State2D(;X,A,B,cells,types)

    if p.epi.loop == 1
        N = num_cells(s)
        for i = 1:N-1
            add_edge!(s.apicalcons, i, i+1, get_type(s,i).apical_junction_init)
            add_edge!(s.basalcons, i, i+1)
        end
        add_edge!(s.apicalcons, N, 1, get_type(s,N).apical_junction_init)
        add_edge!(s.basalcons, N, 1)
        
    else 
        if p.epi.basal_curvature == 0.0
            init_apical_junctions(s, p.epi.init_apical_junction_dist + eps(10.))
            init_basal_junctions(s, p.epi.init_basal_junction_dist + eps(10.))
        else 
            N = num_cells(s)
            for i = 1:N-1
                add_edge!(s.apicalcons, i, i+1, get_type(s,i).apical_junction_init)
                add_edge!(s.basalcons, i, i+1)
            end
        end
    end
    init_cell_events!(s, p)

    return s
end


# centered range
function crange(center::Float64, step::Float64; length::Int64)
    interval_length = step*(length-1)
    return range(center - 0.5*interval_length, center + 0.5*interval_length, length=length)
end

function init_positions(p::StaticParameters)
    epi = p.epi
    N = epi.N_init
    

    rect_domain = epi.init_zone
    left, right = get_xlims(rect_domain)
    bottom, top = get_ylims(rect_domain)

    cells_rect = deepcopy(rect_domain)
    cells_rect.size = cells_rect.size .* (1., 2. / 3.)
    cells_bottom, cells_top = get_ylims(cells_rect)

    X = zeros(2, N)
    A = zeros(2, N)
    B = zeros(2, N)

    center_x = 0.5*(left+right)

    dx = (right-left)/(N+1)
    for (i, x) in enumerate(crange(center_x, dx, length=N))
        X[1,i] = x + rand(Uniform(-dx,dx))
        X[2,i] = rand(Uniform(cells_bottom, cells_top))
    end

    da = epi.init_apical_junction_dist
    for (i, a) in enumerate(crange(center_x, da, length=N))
        A[1,i] = a
        A[2,i] = top
    end

    db = epi.init_basal_junction_dist
    for (i, b) in enumerate(crange(center_x, db, length=N))
        B[1,i] = b
        B[2,i] = bottom
    end

    Rc = 1 / epi.basal_curvature 
    Bc = 1 / epi.basal_curvature_elliSEMTor_factor * Rc

    if epi.basal_curvature  != 0.0

        L = (maximum(B[1,:]) - minimum(B[1,:])) * (1 + 1/N)

        function transform!(x)
            h = x[2,:] .- cells_bottom
            x .-= (0.0, Rc)
            x[1,:] .*= 2*pi/L
            @. x[2,:] = (h + abs(Rc)) * cos(x[1,:] + epi.init_curve_start)
            @. x[1,:] = (h + abs(Bc)) * sin(x[1,:] + epi.init_curve_start)
            x .+= (0.0, Rc)
        end

        transform!(X)
        transform!(A)
        transform!(B)
    end 

    #state.X = HArray{Float64,2}(X)
    #state.A = HArray{Float64,2}(A)
    #state.B = HArray{Float64,2}(B)
    return X, A, B
end


function generate_init_cell_types(p::StaticParameters)
    @unpack N_init, init_distr, init_method = p.epi

    v = split.(split(init_distr, ","), ":")

    types_ordered = [string(strip(type_prob[1])) for type_prob in v]
    probs = [abs( parse(Float64,type_prob[2]) ) for type_prob in v]
    probs_sum = sum(probs)
    # make sure that the probabilities sum up to 1.
    probs ./= probs_sum

    # shuffle the indices, if needed
    indices = (init_method == random) ? shuffle(1:N_init) : collect(1:N_init)

    # generate the types
    cell_type_names_array = Array{String}(undef,N_init)

    cur_type_idx = 1
    p = 0.5 / N_init
    for i in indices
        if p > probs[cur_type_idx]
            p = 0.5 / N_init
            cur_type_idx += 1
        end

        cell_type_names_array[i] = types_ordered[min(length(types_ordered),cur_type_idx)]
        p += 1.0 / N_init
    end
    return cell_type_names_array
end

function init_apical_junctions(s::State, distance )
    N = num_cells(s)
    for j = 1:N
        for i = 1:j-1
            if norm(s.A[:,i] - s.A[:,j]) <= 1.1*distance
                add_edge!(s.apicalcons, i, j, get_type(s,i).apical_junction_init)
            end
        end
    end
end

function init_basal_junctions(s::State, distance )
    N = num_cells(s)
    for j = 1:N
        for i = 1:j-1
            if norm(s.A[:,i] - s.A[:,j]) <= 1.1*distance
                add_edge!(s.basalcons, i, j)
            end
        end
    end
end

function init_cell_events!(s::State, p::StaticParameters)
    for i in 1:num_cells(s)
        for event in values(get_type(s,i).cell_events)
            init_cell_event!(s, p, i, event)
        end
    end
end






#=
 3D
=#
function init_state_3d(p::StaticParameters, rep = 1)
    Random.seed!(p.stat.random_seed + rep - 1)
    X, A, B = init_positions_3d(p)

    init_cell_type_names = generate_init_cell_types(p)
    all_cell_types = [cts.prototype.name for cts in p.cell_types]
    N = p.epi.N_init

    types = EpiCellType[]
    cells = StructVector(EpiCell[])

    for idx in 1:N
        if !(init_cell_type_names[idx] in all_cell_types)
            @error """Unknown type $(init_cell_type_names[idx]) in 'init_distr'. Could it be that 'init_distr' is not correct? \n
            ('init_distr' is $(p.epi.init_distr); known cell types are $(all_cell_types).)"""
        end
        new_type = generate_cell_type(p, init_cell_type_names[idx])
        push!(types, new_type)

        type_idx = length(types)
        epi_cell = EpiCell(new_type; label = idx, type_idx=type_idx)
        init_cell!(epi_cell, new_type)
        push!(cells, epi_cell)
    end

    s = State3D(;X,A,B,cells,types)

    #init_apical_junctions(s, p.epi.init_apical_junction_dist + eps(10.))
    #init_basal_junctions(s, p.epi.init_basal_junction_dist + eps(10.))

    init_cell_events!(s, p)
    return s
end


function init_positions_3d(p::StaticParameters)
    epi = p.epi
    N = epi.N_init

    rect_domain = epi.init_zone
    left, right = get_xlims(rect_domain)
    bottom, top = get_ylims(rect_domain)

    X = zeros(2, N)
    A = zeros(2, N)
    B = zeros(2, N)

    center_x = 0.5*(left+right)

    for i in 1:N
        X[1,i] = rand(Uniform(left,right))
        X[2,i] = rand(Uniform(bottom, top))
        X[3,i] = rand(Uniform(left,right))

        A[1,i] = X[1,i]
        A[2,i] = top
        A[3,i] = X[3,i]

        B[1,i] = X[1,i]
        B[2,i] = bottom
        B[3,i] = X[3,i]
    end

    return X, A, B
end
