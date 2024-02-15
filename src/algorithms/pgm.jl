function sub_time_steps!(s::State, alg::PGM)

    sc = s.cache
    # collision detection
    dt = s.params.sim.dt
    dtf = dt / alg.n_substeps
    μi = 1. / s.params.epi.mu

    for k = 1:alg.n_substeps

        # add forces
        compute_forces!(s)
        s.X .+= μi * dtf * sc.forceX
        s.A .+= μi * dtf * sc.forceA
        s.B .+= μi * dtf * sc.forceB

        # fix constraints
        update_lists!(sc.verletlist, s.X)
        for (i,j) in VerletPairs(sc.verletlist)
            Rᵢⱼ = s[i].R_hard + s[j].R_hard
            sqd = sqdist(s.X[:,i], s.X[:,j])
            if sqd < Rᵢⱼ^2
                # handle overlapping constraints
                d = sqrt(sqd)
                cᵢⱼ = d - Rᵢⱼ  # should be ≧ 0
                Dcᵢ = (s.X[:,i] - s.X[:,j]) / d  #points from j to i
                denom = 2. + alg.alpha
                s.X[:,i] .-= Dcᵢ * cᵢⱼ / denom
                s.X[:,j] .+= Dcᵢ * cᵢⱼ / denom

                # force correction
                s.X[:,i] .+= 0.5 * μi * dtf * (sc.forceX[:,j] + sc.forceX[:,i])
                s.X[:,j] .+= 0.5 * μi * dtf * (sc.forceX[:,j] + sc.forceX[:,i])
            end
        end

        for (i,j) in s.basalcons.edges
            cᵢⱼ = s.B[1,j] - s.B[1,i]

            if cᵢⱼ < 0
                # we only work on the x-coordinate
                Dcᵢ = -1
                Dcⱼ = +1
                denom = 2. + alg.alpha
                s.B[1,i] += Dcᵢ * (-cᵢⱼ) / denom
                s.B[1,j] += Dcⱼ * (-cᵢⱼ) / denom

                # force correction
                s.B[1,i] += 0.5 * μi * dtf * (sc.forceB[1,j] + sc.forceB[1,i])
                s.B[1,j] += 0.5 * μi * dtf * (sc.forceB[1,j] + sc.forceB[1,i])
            end
        end

        for (i,j) in s.basalcons.edges
            v = s.B[1,i] - s.B[1,j]
            max_dist = 0.5*s.cells[i].max_basal_junction_dist +
                        0.5*s.cells[j].max_basal_junction_dist
            cᵢⱼ = max_dist - abs(v)

            if cᵢⱼ < 0
                # we only work on the x-coordinate
                Dcᵢ = -sign(v)
                Dcⱼ = -Dcᵢ
                denom = 2. + alg.alpha
                s.B[1,i] += Dcᵢ * (-cᵢⱼ) / denom
                s.B[1,j] += Dcⱼ * (-cᵢⱼ) / denom

                # force correction
                s.B[1,i] += 0.5 * μi * dtf * (sc.forceB[1,j] + sc.forceB[1,i])
                s.B[1,j] += 0.5 * μi * dtf * (sc.forceB[1,j] + sc.forceB[1,i])
            end
        end

        s.B[2,s.cache.b_constr] .= s.cache.b_constr_value
    end
end
