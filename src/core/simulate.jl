function simulate(p::StaticParameters, rep::Int64 = 1; obs = missing, inplace=false, verbose=nothing, 
    obs_callback=missing, pbd=true)
    s = init_state(p, rep)

    if !isnothing(verbose)
        println("Finished initialisation.")
    end

    return simulate(s, p, rep; obs, inplace, verbose, obs_callback, pbd)
end

function simulate_ensemble(p, n::Int64 = p.stat.num_rep; rep_start = p.stat.random_seed, obs = missing, inplace=false, verbose=nothing, obs_callback=missing, pbd=true)
    ens = Vector{Vector{State{2}}}(undef, n)
    Threads.@threads for i in 1:n
        ens[i] = simulate(p, i; obs, inplace, verbose, obs_callback, pbd)
    end
    return ens
end

function simulate(s::State{Dim}, p::StaticParameters, rep::Int64 = 1; 
    obs = missing, inplace=false, verbose=nothing, obs_callback=missing, pbd=true) where {Dim}

    @unpack epi, sim, alg = p
    @unpack t_end = sim
    @unpack dt = alg
    @unpack init_zone = epi

    c = s.cells

    # now we initialise all the variables we need
    sqrtdt = sqrt(dt)
    dtf = dt/alg.n_substeps
    μi = 1.0 / epi.mu

    @assert sim.dt >= dt
    skp_frames = round(Int64, sim.dt/dt)

    no_curvature = (epi.basal_curvature == 0.)
    Rb = 1. / epi.basal_curvature
    Bb = epi.basal_curvature_elliSEMTor_factor
    y0 = init_zone.center[2] - 0.5*init_zone.size[2]
    basal_info = (;no_curvature,Rb,y0,Rc=@SVector[init_zone.center[1], y0 + Rb ])

    rect_max = deepcopy(epi.init_zone)
    rect_max.size = rect_max.size .* (4., 4.)
    R_max = 2 * maximum( s.cells.R_soft )
    #vl = VerletList{Dim}(rect_max.center - 0.5*rect_max.size, rect_max.center + 0.5*rect_max.size, R_max)

    N = num_cells(s)
    forceX = hzeros(Val(Dim),N)
    forceA = hzeros(Val(Dim),N)
    forceB = hzeros(Val(Dim),N)

    if !pbd
        Xn = copy(s.X)
        An = copy(s.A)
        Bn = copy(s.B)

        Xl = copy(s.X)
        Al = copy(s.A)
        Bl = copy(s.B)

        λ_xx = zeros(N, N)
        λ_bmax = zeros(N)
        λ_bmin = zeros(N)

        dp = p.daha
    end 
    
    da = MVector{2,Float64}(zeros(2));
    dx = MVector{2,Float64}(zeros(2));
    db = MVector{2,Float64}(zeros(2));
    xixj = MVector{2,Float64}(zeros(2));

    vl = nothing
    cache = (;forceX,forceA,forceB,da,dx,db,xixj,vl)

    # let the fun begin

    n_steps = ceil(Int64, t_end / dt)
    n_save_steps = ceil(Int64, n_steps / skp_frames )

    if !inplace
        states = Vector{State{Dim}}(undef, n_save_steps)
    end

    #GLMakie.activate!()
    #f = Figure()
    #s_nd = plot_state!(f[1,1], s, p)


    if !isnothing(verbose)
        println("Start first simulation step.")
    end

    N_last = N
    idx_states = 1
    for l = 1:n_steps
        #s_nd[] = s
        #display(f)

        


        s.t = s.t + dt

        s.cells.age .+= dt

        Rb += dt * p.epi.curve_growth
        Rc = @SVector[init_zone.center[1], y0 + Rb ]

        # first, we take care of the events
        @inbounds for i in 1:N
            type = get_type(s,i)
            for event in type.cell_events
                happend = update!(s, p, i, event)

                if !ismissing(obs) && happend
                    obs_callback(s, p, i, event, obs)
                end
            end
            for event in type.special_cell_events
                happend = update!(s, p, i, event)

                if !ismissing(obs) && happend
                    obs_callback(s, p, i, event, obs)
                end
            end
        end

        c.on_basal_layer[:] .= false
        @inbounds for edge in s.basalcons.edges
            c.on_basal_layer[edge[1]] = true
            c.on_basal_layer[edge[2]] = true
        end

        # make sure the cache is up to date
        N = num_cells(s)
        if N != N_last
            # we have to resize the force vectors
            forceX    = hzeros(Val(Dim),N)
            forceA    = hzeros(Val(Dim),N)
            forceB    = hzeros(Val(Dim),N)
        
            if !pbd
                Xn = copy(s.X)
                An = copy(s.A)
                Bn = copy(s.B)

                Xl = copy(s.X)
                Al = copy(s.A)
                Bl = copy(s.B)
                        
                λ_xx = zeros(N, N)
                λ_bmax = zeros(N)
                λ_bmin = zeros(N)
            end 
            
            cache     = (;forceX,forceA,forceB,da,dx,db,xixj,vl)
        end
        N_last = N


        # new, we update the springs
        time_step_springs!(s, p, basal_info)

        # add noise
        @inbounds @simd for i = 1:N
            f = sqrtdt * sqrt(2.0*s.cells.diffusion[i])
            s.X[:,i] .+= f * @SVector randn(2)
        end

        #
        #adapt_bounds!(vl, s)

        # PBD
        if pbd || s.t < 0.5
            # collision detection

            for k = 1:p.alg.n_substeps
                # add forces
                compute_forces!(s,cache)

                @inbounds for i in 1:N
                    s.X[:,i] .+= μi * dtf * cache.forceX[:,i]
                    s.A[:,i] .+= μi * dtf * cache.forceA[:,i]
                    if !c.is_running[i]
                        s.B[:, i] .+= μi/c.basal_damping_ratio[i] * dtf * cache.forceB[:, i]
                    end
                end


                # fix constraints
                #update_lists!(vl, s.X)
                #@inbounds for (i,j) in VerletPairs(vl)
                @inbounds for i = 1:N
                    for j = 1:i-1
                        Rᵢⱼ = c.R_hard[i] + c.R_hard[j]
                        sqd = sqdist(s.X[:,i], s.X[:,j])
                        if sqd < Rᵢⱼ^2
                            # handle overlapping constraints
                            d = sqrt(sqd)
                            cᵢⱼ = d - Rᵢⱼ  # should be ≧ 0
                            Dcᵢ = (s.X[:,i] - s.X[:,j])  #points from j to i
                            s.X[:,i] -= cᵢⱼ/(2*d) * Dcᵢ
                            s.X[:,j] += cᵢⱼ/(2*d) * Dcᵢ
                        end
                    end
                end

                if no_curvature
                    @inbounds for (i,j) in s.basalcons.edges
                        cᵢⱼ = s.B[1,j] - s.B[1,i]

                        if cᵢⱼ < 0
                            # we only work on the x-coordinate
                            Dcᵢ = -1
                            Dcⱼ = +1
                            s.B[1,i] += Dcᵢ * (-cᵢⱼ) / 2
                            s.B[1,j] += Dcⱼ * (-cᵢⱼ) / 2
                        end
                    end
                else
                    # do nothing
                end

                @inbounds for (i,j) in s.basalcons.edges
                    v = dist(s.B[:,i], s.B[:,j])
                    max_dist = 0.5*c.max_basal_junction_dist[i] +
                                0.5*c.max_basal_junction_dist[j]
                    cᵢⱼ = max_dist - v

                    if cᵢⱼ < 0
                        # we only work on the x-coordinate
                        Dcᵢ = (s.B[:,j] - s.B[:,i]) ./ v
                        Dcⱼ = -Dcᵢ
                        s.B[:,i] .+= Dcᵢ * (-cᵢⱼ) / 2
                        s.B[:,j] .+= Dcⱼ * (-cᵢⱼ) / 2
                    end
                end

                if no_curvature
                    @inbounds for i in 1:N
                        if c.on_basal_layer[i]
                            s.B[2,i] = y0
                        end
                    end
                else
                    @inbounds for i in 1:N
                        if c.on_basal_layer[i]
                            # for projection onto a circle 
                            s.B[:,i] .= abs(Rb) * normalize( (Bb, 1.0) .* (s.B[:,i] .- Rc) ) ./ (Bb, 1.0) + Rc
                        end
                    end
                end
            end
        else  #DAHA 
            
            residual = Inf64 
            max_iter = 1e6
            daha_steps = 0

            Xn .= s.X 
            An .= s.A 
            Bn .= s.B 
            Xl .= s.X 
            Al .= s.A 
            Bl .= s.B 

            λ_xx .= 0.0
            λ_bmin .= 0.0
            λ_bmax .= 0.0
                
            dtc = 1+dp.dt*dp.c/2

            while residual > p.daha.eps && daha_steps < max_iter 

                compute_forces!(s,cache)

                # second order time-step and damping 
                @inbounds for i in 1:N
                    Xn[:,i] .= 1/dtc * ( 
                                2*s.X[:,i] - (2-dtc) * Xl[:,i] 
                                + dp.alpha * dp.dt^2 * cache.forceX[:,i] )
                               
                    An[:,i] .=  1/dtc * ( 2*s.A[:,i] - (2-dtc) * Al[:,i] 
                                + dp.alpha * dp.dt^2 * cache.forceA[:,i] )
                               
                    if !c.is_running[i]
                        Bn[:, i] .= 1/dtc * ( 2*s.B[:,i] - (2-dtc) * Bl[:,i] 
                                        + dp.rho * dp.alpha*dp.dt^2 * cache.forceB[:,i] )
                    else 
                        Bn[:, i] .= s.B[:,i]
                    end
                end
                

                # fix constraints
                #update_lists!(vl, s.X)
                #@inbounds for (i,j) in VerletPairs(vl)
                @inbounds for i = 1:N
                    for j = 1:i-1
                        Rᵢⱼ = c.R_hard[i] + c.R_hard[j]
                        sqd = sqdist(s.X[:,i], s.X[:,j])
                        if sqd < Rᵢⱼ^2
                            # handle overlapping constraints
                            d = sqrt(sqd)
                            cᵢⱼ = -(d^2 - Rᵢⱼ^2)  # should be ≤ 0
                            Dcᵢ = -2*(s.X[:,i] - s.X[:,j])  #points from j to i
                            factor = (dp.alpha*dp.dt^2/dtc * λ_xx[i,j]  
                            + dp.gamma*dp.dt^2/dtc * λ_xx[i,j] * cᵢⱼ )
                            Xn[:,i] .-= factor * Dcᵢ
                            Xn[:,j] .+= factor * Dcᵢ

                            λ_xx[i,j] = max(0.0, λ_xx[i,j] + dp.beta * dp.dt * cᵢⱼ)
                        end
                    end
                end

                
                @inbounds for (i,j) in s.basalcons.edges
                    cᵢⱼ = -(s.B[1,j] - s.B[1,i])

                    if cᵢⱼ < 0
                        # we only work on the x-coordinate
                        Dcᵢ = +1
                        Dcⱼ = -1
                        Bn[1,i] -= ( dp.rho*dp.alpha*dp.dt^2/dtc * λ_bmin[i]  
                                    + dp.gamma*dp.dt^2/dtc * λ_bmin[i] * cᵢⱼ ) * Dcᵢ
                        Bn[1,j] -= ( dp.rho*dp.alpha*dp.dt^2/dtc * λ_bmin[i]  
                                    + dp.gamma*dp.dt^2/dtc * λ_bmin[i] * cᵢⱼ ) * Dcⱼ

                        λ_bmin[i] = max(0.0, λ_bmin[i] + dp.beta * dp.dt * cᵢⱼ)
                    end
                end

                
                @inbounds for (i,j) in s.basalcons.edges
                    v = dist(s.B[:,i], s.B[:,j])
                    max_dist = 0.5*c.max_basal_junction_dist[i] +
                                0.5*c.max_basal_junction_dist[j]
                    cᵢⱼ = -(max_dist - v)

                    if cᵢⱼ > 0
                        # we only work on the x-coordinate
                        Dcᵢ = -(s.B[:,j] - s.B[:,i])
                        Dcⱼ = -Dcᵢ
                        
                        Bn[:,i] .-= ( dp.rho*dp.alpha*dp.dt^2/dtc * λ_bmax[i]  
                                      + dp.gamma*dp.dt^2/dtc * λ_bmax[i] * cᵢⱼ ) * Dcᵢ
                        Bn[:,j] .-= ( dp.rho*dp.alpha*dp.dt^2/dtc * λ_bmax[i]  
                                      + dp.gamma*dp.dt^2/dtc * λ_bmax[i] * cᵢⱼ ) * Dcⱼ

                        λ_bmax[i] = max(0.0, λ_bmax[i] + dp.beta * dp.dt * cᵢⱼ)
                    end
                end

                
                if no_curvature
                    @inbounds for i in 1:N
                        if c.on_basal_layer[i]
                            Bn[2,i] = y0
                        end
                    end
                else
                    @inbounds for i in 1:N
                        if c.on_basal_layer[i]
                            Bn[:,i] .= abs(Rb) * normalize(Bn[:,i] .- Rc) + Rc
                        end
                    end
                end

                residual = 0.0 
                rel_norm = 0.0

                @simd for i in 1:N 
                    residual += (Xn[1,i] - s.X[1,i])^2 + (Xn[2,i] - s.X[2,i])^2
                    residual += (An[1,i] - s.A[1,i])^2 + (An[2,i] - s.A[2,i])^2
                    residual += (Bn[1,i] - s.B[1,i])^2 + (Bn[2,i] - s.B[2,i])^2
                    rel_norm += s.X[1,i]^2 + s.X[2,i]^2 + s.A[1,i]^2 + s.A[2,i]^2 + s.B[1,i]^2 + s.B[2,i]^2
                end

                residual = sqrt(residual) / ( dp.dt * sqrt(rel_norm) )
                daha_steps += 1

                Xl .= s.X 
                Al .= s.A 
                Bl .= s.B 
                s.X .= Xn 
                s.A .= An
                s.B .= Bn  
            end 

            if daha_steps == max_iter 
                @error "DAHA did not converge."
                return missing
            end

            
            if !ismissing(obs)
                obs_callback(s, p, 0, (;type=:daha, t=s.t, steps = daha_steps), obs)
            end
        end

        if !inplace && mod(l,skp_frames) == 0

            if !isnothing(verbose)
                verbose(l, s, p, cache)
            end

            states[idx_states] = quickcopy(s)
            idx_states += 1
        end
    end



    if !inplace
        resize!(states, idx_states-1)
        return states
    else
        s
    end
end


function time_step_springs!(s::State, p::StaticParameters, basal_info)
    N = num_cells(s)
    dt = p.alg.dt
    c = s.cells

    #TODO replace running_zone with the correct condition for curved tissues.

    # check which cells are running
    @inbounds for i in 1:N
        c.is_below[i] = is_below_basal_layer(s, p, i, basal_info)
        rm = s.cells.running_mode[i]
        c.is_running[i] = !c.on_basal_layer[i] && s.B[2,i] > p.epi.running_zone[1] && (rm >= 3 ||
                            (c.is_below[i] && rm >= 1))
    end

    @inbounds for i in 1:N
        cell = s[i]
        # update cytoskeletons, 'drl' = desired_rest_length
        if cell.apical_cytos_strain == 0.
            apical_drl = max( 0., dist(s.X[:,i], s.A[:,i]) - cell.R_soft)
        elseif cell.apical_cytos_strain < 0. && cell.apical_cytos_strain >= -1 # movement towards apical layer
            apical_drl = (cell.apical_cytos_strain + 1) * max( 0., dist(s.X[:,i], s.A[:,i]) - cell.R_soft)
        elseif cell.apical_cytos_strain > 0  # movement towards basal layer
            apical_drl = cell.apical_cytos_strain * max( 0., dist(s.A[:,i], s.B[:,i]) - 2. * cell.R_soft) + (cell.apical_cytos_strain - 1) * max( 0., dist(s.X[:,i], s.A[:,i]) - cell.R_soft)
        else
            @error "apical_cytos_strain ∉ {-1,0,1} is not implemented."
        end

        # explicit solution
        cell.apical_cytos_rest_length = exp(-dt * cell.k_cytoskeleton ) *
                (cell.apical_cytos_rest_length - apical_drl) + apical_drl

        #cell.apical_cytos_rest_length = max(cell.apical_cytos_rest_length, cell.max_apical_cytos_rest_length)

        if cell.running_mode == 0 || c.on_basal_layer[i]
            if cell.basal_cytos_strain == 0.
                basal_drl = max( 0., dist(s.X[:,i], s.B[:,i]) - cell.R_soft)
            elseif cell.basal_cytos_strain > 0 && cell.basal_cytos_strain <= 1 # movement towards apical layer
                basal_drl = cell.basal_cytos_strain * max( 0., dist(s.A[:,i], s.B[:,i]) - 2. * cell.R_soft) + (1-cell.basal_cytos_strain) * max( 0., dist(s.X[:,i], s.B[:,i]) - cell.R_soft)
            elseif cell.basal_cytos_strain > 1 
                basal_drl = (cell.basal_cytos_strain-1) * max( 0., dist(s.X[:,i], s.B[:,i]))
            elseif cell.basal_cytos_strain < 0
                basal_drl = (1+cell.basal_cytos_strain) * max( 0., dist(s.X[:,i], s.B[:,i]) - cell.R_soft)
            else
                @error "basal_cytos_strain ∉ {-1,0,1} is not implemented."
            end

            # explicit solution
            cell.basal_cytos_rest_length = exp(-dt * cell.k_cytoskeleton ) *
                (cell.basal_cytos_rest_length - basal_drl) + basal_drl


            # enforce maximal length:
            too_long = cell.basal_cytos_rest_length + cell.apical_cytos_rest_length - p.epi.max_cytoskeleton_length

            if too_long > 0
                cell.basal_cytos_rest_length *=  ( p.epi.max_cytoskeleton_length - too_long ) / p.epi.max_cytoskeleton_length
                cell.apical_cytos_rest_length *= ( p.epi.max_cytoskeleton_length - too_long ) / p.epi.max_cytoskeleton_length
            end

        else
            # running phase
            if (cell.running_mode >= 2 && !c.is_below[i])
                cell.basal_cytos_rest_length += dt*cell.running_speed
            else
                cell.basal_cytos_rest_length =  exp(-dt * cell.k_cytoskeleton ) *
                        (cell.basal_cytos_rest_length)
            end

            if c.is_running[i]
                s.B[:,i] .+= dt*( cell.running_speed * normalize(s.B[:,i] - s.X[:,i]) )
            end
        end

        #cell.basal_cytos_rest_length = max(cell.basal_cytos_rest_length, cell.max_basal_cytos_rest_length)
    end

    ac = s.apicalcons
    for ℓ in eachindex(ac.edges)
        (u,v) = ac.edges[ℓ]

        junction_drl = 0.

        if s.cells.apical_junction_strain[u] == 0. && s.cells.apical_junction_strain[v] == 0.
            junction_drl = max( 0., dist(s.A[:,u], s.A[:,v]))
        elseif s.cells.apical_junction_strain[u] == -1. || s.cells.apical_junction_strain[v] == -1.
            junction_drl = 0.
        else
            @error "apical_junction_strain ∉ {-1,0} is not implemented."
        end

        k_avg = 0.5*s.cells.k_apical_junction[u] + 0.5*s.cells.k_apical_junction[v]
        # explicit solution
        ac.data[ℓ] = exp(-dt * k_avg ) * (ac.data[ℓ] - junction_drl) + junction_drl
    end
end
