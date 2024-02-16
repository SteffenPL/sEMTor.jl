function pbd!(s::State, p::StaticParameters, basal_info, cache)
    (;no_curvature, Rb, y0, Rc) = basal_info
    (;vl) = cache

    c = s.cells
    # collision detection
    dtf = p.sim.dt / p.alg.n_substeps
    μi = 1. / p.epi.mu

    N = size(s.X,2)::Int64

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
        update_lists!(vl, s.X)
        @inbounds for (i,j) in VerletPairs(vl)
            Rᵢⱼ = c.R_hard[i] + c.R_hard[j]
            sqd = sqdist(s.X[:,i], s.X[:,j])
            if sqd < Rᵢⱼ^2
                # handle overlapping constraints
                d = sqrt(sqd)
                cᵢⱼ = d - Rᵢⱼ  # should be ≧ 0
                Dcᵢ = (s.X[:,i] - s.X[:,j]) / d  #points from j to i
                denom = 2. + p.alg.alpha
                s.X[:,i] .-= Dcᵢ * cᵢⱼ / denom
                s.X[:,j] .+= Dcᵢ * cᵢⱼ / denom
            end
        end

        @inbounds for (i,j) in s.basalcons.edges
            cᵢⱼ = s.B[1,j] - s.B[1,i]

            if cᵢⱼ < 0
                # we only work on the x-coordinate
                Dcᵢ = -1
                Dcⱼ = +1
                denom = 2. + p.alg.alpha
                s.B[1,i] += Dcᵢ * (-cᵢⱼ) / denom
                s.B[1,j] += Dcⱼ * (-cᵢⱼ) / denom
            end
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
                denom = 2. + p.alg.alpha
                s.B[:,i] .+= Dcᵢ * (-cᵢⱼ) / denom
                s.B[:,j] .+= Dcⱼ * (-cᵢⱼ) / denom
            end
        end

        if no_curvature
            s.B[2,c.on_basal_layer] .= y0
        else
            @inbounds for i in 1:N
                if c.on_basal_layer[i]
                    s.B[:,i] .= abs(Rb) * normalize(s.B[:,i] .- Rc) + Rc
                end
            end
        end
    end
end
