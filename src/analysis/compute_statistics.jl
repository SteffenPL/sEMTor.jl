function project_onto_line(L, P)
    v_L = L[:,2] - L[:,1]
    v_P = P   - L[:,1]
    d_L = dot(v_L, v_L)
    if d_L == 0.
        return L[:,1]
    end
    f = max(0., min(1., dot(v_L,v_P) / dot(v_L, v_L)))
    return L[:,1] + f * v_L
end

function project_onto_linestrip(L, P)
    dist = Inf64
    Pr = [NaN, NaN]
    for i in 2:size(L,2)
        Pr_i = project_onto_line(L[:,i-1:i], P)
        dist_i = norm(Pr_i - P)
        if dist_i < dist
            Pr = Pr_i
            dist = dist_i
        end
    end
    return Pr
end

function project_onto_apical_layer(s, i)
	has_apical_ind = falses(num_cells(s))
	for (u,v) in s.apicalcons.edges
		has_apical_ind[u] = true
		has_apical_ind[v] = true
	end
	return project_onto_linestrip(s.A[:,has_apical_ind], s.X[:,i])
end

function project_onto_basal_layer(s, i)
	has_basal_ind = falses(num_cells(s))
	for (u,v) in s.basalcons.edges
		has_basal_ind[u] = true
		has_basal_ind[v] = true
	end
	return project_onto_linestrip(s.B[:,has_basal_ind], s.X[:,i])
end

function compute_statistics(state_reps::Vector{Vector{State2D}}, p::StaticParameters, stats_dt=p.sim.dt; control_name = "control")

	if( length(state_reps) == 0 || length(state_reps[1]) == 0)
		return missing
	end

	n_reps = length(state_reps)
	ct_names = vcat("tissue", cell_type_names(p))

	# make sure ct_names follows "tissue":
	if control_name in ct_names
		control_type_idx = findfirst(isequal(control_name), ct_names)
		# swap
		ct_names[control_type_idx] = ct_names[2]
		ct_names[2] = control_name
	end
	n_types = length(ct_names)

	no_curvature = (p.epi.basal_curvature == 0.)
    Rb = 1. / p.epi.basal_curvature
    y0 = p.epi.init_zone.center[2] - 0.5*p.epi.init_zone.size[2]
    Rc = @SVector[p.epi.init_zone.center[1], y0 + Rb ]
    basal_info = (;no_curvature,Rb,y0,Rc)


	dt = p.sim.dt
	detail = floor(Int64, stats_dt / dt)

    columns = [:time, :rep, 
               :apical_width, :nuclei_width, :basal_width, :apical_basal_width_ratio,
               :basal_apical_height, :basal_nuclei_height, :nuclei_apical_height,
			   :apical_height, :nuclei_height, :basal_height,
               :apical_straightness, :cell_count, :center_count,
               :basal_extrusion_count, :basal_extrusion_ratio,
               :apical_extrusion_count, :apical_extrusion_ratio,
               :below_control_count, :below_control_ratio]

	n_steps = length(1:detail:length(state_reps[1]))
	n_rows = n_reps * n_steps

   	stats = Dict(
   				ct_name => DataFrame(Dict(col=>zeros(n_rows) for col in columns))
   				for ct_name in ct_names
   				)
   	div_stats = Dict(
   				ct_name => DataFrame()
   				for ct_name in ct_names
   				)

	for rep in 1:n_reps
		states = state_reps[rep][1:detail:end]

		rows_start = (rep-1)*n_steps+1
		for type in ct_names
			stats[type].rep[rows_start:rows_start+n_steps-1] .= rep
		end

    	M = length(states)

		if M == 0
			continue
		end

		N = num_cells(states[1])
		type_ind = falses(N)
		control_ind = falses(N)
		has_apical_ind = falses(N)
		has_basal_ind = falses(N)
	    center_ind = falses(N)

	    for (t_idx, s) in enumerate(states)

			row = t_idx + rows_start-1
			x1min =  minimum(s.X[1,:])
			x1max = maximum(s.X[1,:])
			N = num_cells(s)
			if N > length(center_ind)
				resize!(center_ind, N)
				resize!(control_ind, N)
				resize!(has_apical_ind, N)
				resize!(has_basal_ind, N)
				resize!(type_ind, N)
			end

			for k in 1:N
				control_ind[k] = get_type(s,k).name == control_name
			end
			x2min_control = minimum(s.X[2,control_ind])

			has_apical_ind .= false
			for (u,v) in s.apicalcons.edges
				has_apical_ind[u] = true
				has_apical_ind[v] = true
			end

			has_basal_ind .= false
			for (u,v) in s.basalcons.edges
				has_basal_ind[u] = true
				has_basal_ind[v] = true
			end

			df_tissue = stats["tissue"]

			for (j, type) in enumerate(ct_names)
				df = stats[type]
				df.time[row] = s.t

				for k in 1:N
					type_ind[k] = (type == "tissue" || get_type(s,k).name == type)
				end

		        df.cell_count[row] = count(type_ind)
				if df.cell_count[row] == 0
					continue
				end

		        df.apical_width[row] = maximum(s.A[1,type_ind]) - minimum(s.A[1,type_ind])
		        df.nuclei_width[row] = maximum(s.X[1,type_ind]) - minimum(s.X[1,type_ind])
		        df.basal_width[row] = maximum(s.B[1,type_ind]) - minimum(s.B[1,type_ind])
				df.apical_basal_width_ratio[row] = df.apical_width[row] / df.basal_width[row]


				aw = df_tissue.nuclei_width[row]  # tissue width
				for k in 1:N
					center_ind[k] = type_ind[k] && (x1min + 0.25*aw <= s.X[1,k] <= 0.25*aw + x1max)
				end

				df.center_count[row] = count(center_ind)

				#=
				df.basal_apical_height[row]  = mean(s.A[2,center_ind] .- s.B[2,center_ind])
				df.basal_nuclei_height[row]  = mean(s.X[2,center_ind] .- s.B[2,center_ind])
				df.nuclei_apical_height[row] = mean(s.A[2,center_ind] .- s.X[2,center_ind])
				=#

				df.apical_height[row]  = mean(s.A[2,type_ind])
				df.nuclei_height[row]  = mean(s.X[2,type_ind])
				df.basal_height[row]  = mean(s.B[2,type_ind])

				df.basal_apical_height[row]  = mean(s.A[2,type_ind] .- s.B[2,type_ind])
				df.basal_nuclei_height[row]  = mean(s.X[2,type_ind] .- s.B[2,type_ind])
				df.nuclei_apical_height[row] = mean(s.A[2,type_ind] .- s.X[2,type_ind])


			    a_length = 0.
				ind_left = (1:N)[type_ind][findmin(s.B[1,type_ind])[2]]
				# this is a quick hack to translate the index w.r.t.
				# to the type back to an index w.r.t. the tissue
				ind_right = (1:N)[type_ind][findmax(s.B[1,type_ind])[2]]
				for (u, v) in s.apicalcons.edges
					if type_ind[u] && type_ind[v]
						a_length +=  SEMTor.dist( s.A[:,u], s.A[:,v])
					end
				end
				df.apical_straightness[row] = SEMTor.dist(s.A[:,ind_left], s.A[:,ind_right]) / a_length

				df.basal_extrusion_count[row] = count( nuclei_is_below_basal_layer(s, p, i, basal_info) for i in 1:N if type_ind[i] )
				df.below_control_count[row] = count(s.X[2,type_ind] .< x2min_control)

				if(s.t == 0.0)
					df.basal_extrusion_count[row] = 0
				end

				apical_extrusion_count = 0
				for k = 1:N
					if type_ind[k]
						if has_apical_ind[k]
							if s.X[2,k] > s.A[2,k]
								apical_extrusion_count += 1
							end
						else
							A_eff = project_onto_linestrip(s.A[:,has_apical_ind], s.X[:,k])
							if s.X[2,k] > A_eff[2]
								apical_extrusion_count += 1
							end
						end
					end
				end
				df.apical_extrusion_count[row] = apical_extrusion_count


				df.basal_extrusion_ratio[row] = df.basal_extrusion_count[row] / df.cell_count[row]
				df.below_control_ratio[row] = df.below_control_count[row] / df.cell_count[row]
				df.apical_extrusion_ratio[row] = df.apical_extrusion_count[row] / df.cell_count[row]
			end 
		end 


		# division statistics need larger dt
		states_full = state_reps[rep]
    	M = length(states_full)

		if M == 0
			continue
		end

		N = num_cells(states_full[1])
		type_ind = falses(N)
		control_ind = falses(N)
		has_apical_ind = falses(N)
		has_basal_ind = falses(N)
	    center_ind = falses(N)

		for (t_idx, s) in enumerate(states_full)

			
			N = num_cells(s)
			if N > length(type_ind)
				resize!(has_apical_ind, N)
				resize!(has_basal_ind, N)
				resize!(type_ind, N)
			end

			for (j, type) in enumerate(ct_names)

				for k in 1:N
					type_ind[k] = (type == "tissue" || get_type(s,k).name == type)
				end
					
				has_apical_ind .= false
				for (u,v) in s.apicalcons.edges
					has_apical_ind[u] = true
					has_apical_ind[v] = true
				end

				has_basal_ind .= false
				for (u,v) in s.basalcons.edges
					has_basal_ind[u] = true
					has_basal_ind[v] = true
				end
				
				# has_apical_ind
				# has_basal_ind
				# type_ind

				df_divisions = div_stats[type]
				if t_idx + 1 <= M
					for k = 1:N
						if type_ind[k] == true
							if s.cells.birth_time[k] < states_full[t_idx+1].cells.birth_time[k]

								Aeff = project_onto_linestrip(s.A[:,has_apical_ind], s.X[:,k])
								Beff = project_onto_linestrip(s.B[:,has_basal_ind], s.X[:,k])
								append!(df_divisions,
									DataFrame(  time=[s.t],
												rep=rep,
												type=type,
												label=[s.cells.label[k]],
												idx=[k],
												X1=[s.X[1,k]],
												X2=[s.X[2,k]],
												division_height=[(s.X[2,k] - Beff[2])/(Aeff[2] - Beff[2])]) )
							end
						end
					end
				end
			end
		end
	end

	stats_mean = Dict()
	for (type, df) in stats
		gb = groupby(df, :time)
		df_mean = combine(gb, valuecols(gb) .=> mean, valuecols(gb) .=> std)
		stats_mean[type] = df_mean
	end

	return (;stats, div_stats, stats_mean)
end
