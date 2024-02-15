function minimal_plot_rect(states::Array{State{Dim}}, p::StaticParameters; margin = 0.2)::Rectangle{Dim} where {Dim}
	plot_rect=deepcopy(p.epi.init_zone)

	for state in states
    	extend!(plot_rect, state.X)
    	extend!(plot_rect, state.A)
    	extend!(plot_rect, state.B)
	end

    plot_rect.size *= 1. + margin

    return plot_rect
end

function minimal_plot_rect(state::State{Dim}, p::StaticParameters; margin = 0.2)::Rectangle{Dim} where {Dim}
	return minimal_plot_rect( [state], p, margin=margin)
end

function get_celltype_colormap(state, p)
    colors = String[]
    colors_transp = Tuple{String,Float64}[]
    type2idx = Dict{String,Int64}()
    i = 1
    for cts in p.cell_types
		type = cts.prototype
        col = type.color
        if col in ("none","")
            col = "green"
        end
        push!(colors, col)
        push!(colors_transp, (col, 0.2))
        type2idx[type.name] = i
        i += 1
    end
    colors, colors_transp, type2idx
end


function plot_state!(f, state, p, rect = minimal_plot_rect(state, p); 
    s_nd = Observable(state), title = p.sim.name * "\n", 
    rep_nd=missing, show_labels=true, show_idxs=false, 
    lw = 0.3, hard_col = false, ms = .2, show_age = true)

    cm, cm_transp, type2idx = get_celltype_colormap(state, p)
    age_alpha = 0.05

    cm_age_hard = [  
        [(1-a)*RGBf(0.9,0.8,0.3) + a*RGBf(0.1,0.6,0.2) 
        for a in LinRange(0,1,10)];
        RGBf(0.05,0.3,0.62);
        RGBf(0.8,0.1,0.1);]
    
    cm_age = [RGBAf(c, age_alpha) for c in cm_age_hard]
        
    nuclei_pos = lift( s -> Point2f.(s.X[1,:],s.X[2,:]), s_nd )
    apical_pos  = lift( s -> Point2f.(s.A[1,:], s.A[2,:]), s_nd)
    basal_pos  = lift( s -> Point2f.(s.B[1,:], s.B[2,:]), s_nd)
    R_hard = lift( s -> 2*s.cells.R_hard, s_nd )
    R_soft = lift( s -> 2*s.cells.R_soft, s_nd )
    cell_age = lift( s -> s.cells.age, s_nd)

    function age_scale(age,div_age) 
        x = age - div_age + 2
        if x > 0 
            return x / 2 
        else 
            return x / (div_age-2)
        end 
    end

    cell_age_rel = lift(s_nd) do  s
        age_scale.(s.cells.age, s.cells.division_time)     
    end

    t = @lift ($s_nd).t
    colors = @lift $s_nd.cells.prototype_idx
    if !ismissing(rep_nd)
        title_str = @lift string(title, "  (rep: ", $rep_nd, ")", "\n t = ", round($t,digits=2), " h")
    else
        title_str = @lift string(title , "t = ", round($t,digits=2), " h")
    end
    apical_nuclei = lift( s_nd ) do s
        [(Point2f(s.X[:,j]),
          Point2f(s.A[:,j]))
          for j in 1:size(s.X,2)]
    end
    apical_cons = lift( s_nd ) do s
        [(Point2f(s.A[:,u]),
          Point2f(s.A[:,v]))
          for (u,v) in s.apicalcons.edges]
    end
    basal_nuclei = lift( s_nd ) do s
        [(Point2f(s.X[:,j]),
          Point2f(s.B[:,j]))
          for j in 1:size(s.X, 2)]
    end
    basal_cons = lift( s_nd ) do s
        [(Point2f(s.B[:,u]),
          Point2f(s.B[:,v]))
          for (u,v) in s.basalcons.edges]
    end


    ax = Axis(f, title = title_str) 
    ax.xrectzoom[] = false
    ax.yrectzoom[] = false

    limits!(ax, SEMTor.get_xlims(rect), SEMTor.get_ylims(rect))


    linesegments!(ax, apical_nuclei, color = :orange, linewidth = lw)
    linesegments!(ax, basal_nuclei, color = :orange, linewidth = lw)
    linesegments!(ax, apical_cons, color=:red, linewidth = lw)
    linesegments!(ax, basal_cons, color=:black, linewidth = lw)

    cell_to_ap_col(c) = c.has_apical_cytos ? :red : :black

    ap_m = lift( s_nd ) do s
        [ any( i in e for e in s.apicalcons.edges) ? :circle : :xcross
          for i in 1:length(s)]
    end

    
    ba_m = lift( s_nd ) do s
        [ any( i in e for e in s.basalcons.edges) ? :circle : :xcross
          for i in 1:length(s)]
    end
    

    scatter!(ax, nuclei_pos,
        markersize = R_soft, markerspace=:data,
        color = cell_age_rel, colormap=cm_age, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
        strokecolor=:black, strokewidth = 0.5)
        
    if show_age
            
        scatter!(ax, nuclei_pos,
        markersize = R_soft, markerspace=:data,
        color = cell_age_rel, colormap=cm_age, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
        strokecolor=:black, strokewidth = 0.5)
        
        for p in LinRange(0.3, 1.0, 10)
            scatter!(ax, nuclei_pos,
            markersize = @lift( p .* $R_soft ), markerspace=:data,
            color = cell_age_rel, colormap=cm_age, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
            strokecolor=(:black), strokewidth = 0.0)
        end
    else 

        scatter!(ax, nuclei_pos,
        markersize = R_soft, markerspace=:data,
        color = colors, colormap=cm_transp, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
        strokecolor=:black, strokewidth = 1.0)
        
    end

    if hard_col
        scatter!(ax, nuclei_pos,
                    markersize = R_hard, markerspace=:data,
                    color = cell_age_rel, colormap = cm_age_hard,
                strokecolor=:black, strokewidth = 0.5)
    else 
        scatter!(ax, nuclei_pos,
                    markersize = R_hard, markerspace=:data,
                    color = colors, colormap = cm,
                strokecolor=:black, strokewidth = 1.0)
    end

    scatter!(ax, apical_pos, marker = ap_m, markersize = ms, markerspace=:data, color = :red)
    scatter!(ax, basal_pos, marker = ba_m, markersize = ms, markerspace=:data, color = :black)


	if show_labels
	    labels = @lift string.($s_nd.cells.label[1:length($nuclei_pos)])
    	#text!(ax, labels, position = nuclei_pos, fontsize = 0.3, align = (:center, :center))
	elseif show_idxs
	    idx = @lift string.(1:length($nuclei_pos))
    	#text!(ax, idx, position = nuclei_pos, fontsize = 0.3, align = (:center, :center))
	end

    ax.aspect = DataAspect()
    ax.xlabel = "[5μm]"
    ax.ylabel = "[5μm]"


    s_nd
end

function update_plot!(obs, state)
    obs[] = state
end

function plot_state(state, p; kwargs...)
    f = Figure()
    plot_state!(f[1,1], state, p; kwargs...)
    return f 
end








function plot_state_black!(f, state, p, rect = minimal_plot_rect(state, p); 
    s_nd = Observable(state), title = p.sim.name * "\n", 
    rep_nd=missing, show_labels=true, show_idxs=false, 
    lw = 0.3, hard_col = false, ms = .2, show_age = true)

    cm, cm_transp, type2idx = get_celltype_colormap(state, p)
    age_alpha = 0.05

    cm_age_hard = [  
        [(1-a)*RGBf(0.9,0.8,0.3) + a*RGBf(0.1,0.6,0.2) 
        for a in LinRange(0,1,10)];
        RGBf(0.05,0.3,0.62);
        RGBf(0.8,0.1,0.1);]
    
    cm_age = [RGBAf(c, age_alpha) for c in cm_age_hard]
        
    nuclei_pos = lift( s -> Point2f.(s.X[1,:],s.X[2,:]), s_nd )
    apical_pos  = lift( s -> Point2f.(s.A[1,:], s.A[2,:]), s_nd)
    basal_pos  = lift( s -> Point2f.(s.B[1,:], s.B[2,:]), s_nd)
    R_hard = lift( s -> 2*s.cells.R_hard, s_nd )
    R_soft = lift( s -> 2*s.cells.R_soft, s_nd )
    cell_age = lift( s -> s.cells.age, s_nd)

    function age_scale(age,div_age) 
        x = age - div_age + 2
        if x > 0 
            return x / 2 
        else 
            return x / (div_age-2)
        end 
    end

    cell_age_rel = lift(s_nd) do  s
        age_scale.(s.cells.age, s.cells.division_time)     
    end

    t = @lift ($s_nd).t
    colors = @lift $s_nd.cells.prototype_idx
    if !ismissing(rep_nd)
        title_str = @lift string(title, "  (rep: ", $rep_nd, ")", "\n t = ", round($t,digits=2), " h")
    else
        title_str = @lift string(title , "t = ", round($t,digits=2), " h")
    end
    apical_nuclei = lift( s_nd ) do s
        [(Point2f(s.X[:,j]),
          Point2f(s.A[:,j]))
          for j in 1:size(s.X,2)]
    end
    apical_cons = lift( s_nd ) do s
        [(Point2f(s.A[:,u]),
          Point2f(s.A[:,v]))
          for (u,v) in s.apicalcons.edges]
    end
    basal_nuclei = lift( s_nd ) do s
        [(Point2f(s.X[:,j]),
          Point2f(s.B[:,j]))
          for j in 1:size(s.X, 2)]
    end
    basal_cons = lift( s_nd ) do s
        [(Point2f(s.B[:,u]),
          Point2f(s.B[:,v]))
          for (u,v) in s.basalcons.edges]
    end


    ax = Axis(f, title = title_str, axis=(xrectzoom = false, yrectzoom = false))

    limits!(ax, SEMTor.get_xlims(rect), SEMTor.get_ylims(rect))


    linesegments!(ax, apical_nuclei, color = :orange, linewidth = lw)
    linesegments!(ax, basal_nuclei, color = :orange, linewidth = lw)
    linesegments!(ax, apical_cons, color=:red, linewidth = lw)
    linesegments!(ax, basal_cons, color=:black, linewidth = lw)

    cell_to_ap_col(c) = c.has_apical_cytos ? :red : :black

    ap_m = lift( s_nd ) do s
        [ any( i in e for e in s.apicalcons.edges) ? :circle : :xcross
          for i in 1:length(s)]
    end

    
    ba_m = lift( s_nd ) do s
        [ any( i in e for e in s.basalcons.edges) ? :circle : :xcross
          for i in 1:length(s)]
    end
    

    scatter!(ax, nuclei_pos,
        markersize = R_soft, markerspace=:data,
        color = cell_age_rel, colormap=cm_age, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
        strokecolor=:black, strokewidth = 0.5)
        
    if show_age
            
        scatter!(ax, nuclei_pos,
        markersize = R_soft, markerspace=:data,
        color = cell_age_rel, colormap=cm_age, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
        strokecolor=:black, strokewidth = 0.5)
        
        for p in LinRange(0.3, 1.0, 10)
            scatter!(ax, nuclei_pos,
            markersize = @lift( p .* $R_soft ), markerspace=:data,
            color = cell_age_rel, colormap=cm_age, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
            strokecolor=(:black), strokewidth = 0.0)
        end
    else 

        scatter!(ax, nuclei_pos,
        markersize = R_soft, markerspace=:data,
        color = colors, colormap=cm_transp, #cm_age, #(:white,0.4), #colormap = cm=cm_age,
        strokecolor=:black, strokewidth = 1.0)
        
    end

    if hard_col
        scatter!(ax, nuclei_pos,
                    markersize = R_hard, markerspace=:data,
                    color = cell_age_rel, colormap = cm_age_hard,
                strokecolor=:black, strokewidth = 0.5)
    else 
        scatter!(ax, nuclei_pos,
                    markersize = R_hard, markerspace=:data,
                    color = colors, colormap = cm,
                strokecolor=:black, strokewidth = 1.0)
    end

    scatter!(ax, apical_pos, marker = ap_m, markersize = ms, markerspace=:data, color = :red)
    scatter!(ax, basal_pos, marker = ba_m, markersize = ms, markerspace=:data, color = :black)


	if show_labels
	    labels = @lift string.($s_nd.cells.label[1:length($nuclei_pos)])
    	text!(ax, labels, position = nuclei_pos, textsize = 0.3, space = :data, align = (:center, :center))
	elseif show_idxs
	    idx = @lift string.(1:length($nuclei_pos))
    	text!(ax, idx, position = nuclei_pos, textsize = 0.3, space = :data, align = (:center, :center))
	end

    ax.aspect = DataAspect()
    ax.xlabel = "[5μm]"
    ax.ylabel = "[5μm]"


    s_nd
end