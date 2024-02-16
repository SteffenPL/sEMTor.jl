function compute_forces!(s::State{Dim}, cache) where {Dim}
    (;forceX, forceA, forceB, da, dx, db, xixj, vl) = cache
    N = num_cells(s)

    c = s.cells

    # set zero
    forceX .= 0.
    forceA .= 0.
    forceB .= 0.

    @inbounds for i=1:N
        da .= 0
        dx .= 0
        db .= 0

        a = s.A[:,i]
        x = s.X[:,i]
        b = s.B[:,i]

        # apical-nuclei spring
        distAX = dist(x, a)
        distBX = dist(x, b)

        ax = x-a
        if c.has_apical_cytos[i] && distAX != 0
            rl = (c.apical_cytos_rest_length[i] + c.R_soft[i])
            f = 2*c.stiffness_nuclei_apical[i] * (distAX - rl) / (distAX*rl^2)

            dx .-= f * ax
            da .+= f * ax
        end
        # basal-nuclei spring
        bx = x-b
        if c.has_basal_cytos[i] && distBX != 0
            rl = (c.basal_cytos_rest_length[i] + c.R_soft[i])
            f = 2*c.stiffness_nuclei_basal[i] * (distBX - rl) / (distBX*rl^2)
            dx .-= f * bx
            db .+= f * bx
        end

        # straightnessR = A - X
        ax_bx = dot(ax,bx)
        if c.has_apical_cytos[i] && c.has_basal_cytos[i] && ax_bx != 0
            f = c.stiffness_straightness[i] / (distAX*distBX)
            dR = f * ( -bx + ax_bx/distAX^2 * ax )
            dS = f * ( -ax + ax_bx/distBX^2 * bx )
            da .-= dR
            dx .+= dR+dS
            db .-= dS
        end


        forceA[:,i] .+= da
        forceX[:,i] .+= dx
        forceB[:,i] .+= (1.0 / c.basal_damping_ratio[i]) .* db
    end

    #update_lists!(vl, s.X)
    #@inbounds for (i,j) in VerletPairs(vl)
    @inbounds for i = 1:N
        for j = 1:i-1
            xixj = s.X[:,i] - s.X[:,j]
            Rᵢⱼ = c.R_soft[i] + c.R_soft[j]
            # soft-spheres
            sqd = sum(x -> x^2, xixj)
            if sqd <= Rᵢⱼ^2
                d = sqrt(sqd)
                f = (c.stiffness_repulsion[i] + c.stiffness_repulsion[j]) /
                    Rᵢⱼ * ( d/Rᵢⱼ - 1. ) / d
                forceX[:,i] .-= f .* xixj
                forceX[:,j] .+= f .* xixj
            end
        end
    end


    apicons = s.apicalcons
    @inbounds for ((u,v),rest_length) in EdgeIter(apicons)
        au = s.A[:,u]
        av = s.A[:,v]
        auav = au - av
        d = sqrt(sum(x->x^2,auav))
        f = 0.25 * 0.5 * (c.stiffness_apical_apical[u] + c.stiffness_apical_apical[v]) *
                    (d - rest_length) / d
        forceA[:,u] .-= f .* auav
        forceA[:,v] .+= f .* auav
    end

    nothing
end


function compute_forces(s::State)
    compute_forces!(s)
    return (s.cache.forceX, s.cache.forceA, s.cache.forceB)
end
