function categorise(t1,t2,t3,t4, n1 = "A", n2 = "B", n3 = "P", n4 = "R", n0="")
    tup(t,n) = (t, isinf(t) ? n0 : n)
    Ts = [tup(t1,n1), tup(t2,n2), tup(t3,n3), tup(t4,n4)]
    sort!(Ts, by = tn -> tn[1])
    res = prod(tn[2] for tn in Ts)
    return isempty(res) ?  "∅" : res
end

function basal_extr_criterium(x, wnd)
	for i in eachindex(x)
		if (    x[i] <= 0.0
			&&	mean(x[i:end]) <= 0.0
			&&	mean(x[i:min(i+wnd-1,length(x))]) <= 0.0 )
			return i
		end
	end
	return -1
end

function analyse_state(states, p, rep, cell_type, obs)

	s = states[end]
    N = num_cells(s)
	x = Vector{Float64}(undef, length(states))
	wnd = ceil(Int64, 1.0 / p.sim.dt)  # should be at least 1 hr below at first entry.

	has_apical_ind = falses(N)
	for (u,v) in s.apicalcons.edges
		has_apical_ind[u] = true
		has_apical_ind[v] = true
	end

	no_curvature = (p.epi.basal_curvature == 0.)
    Rb = 1. / p.epi.basal_curvature
    y0 = p.epi.init_zone.center[2] - 0.5*p.epi.init_zone.size[2]
    Rc = SEMTor.StaticArrays.@SVector[p.epi.init_zone.center[1], y0 + Rb ]
    basal_info = (;no_curvature,Rb,y0,Rc)

	maybe_no_inm = occursin("noINM", p.sim.name) || occursin("no INM", p.sim.name)

    df = DataFrame(type = String[],
			rep = Int64[], label = Int64[], idx = Int64[],
			INM = Bool[],
			A = Float64[], B = Float64[], P = Float64[], R = Float64[],
			running_mode = Int64[],
			apical_extr = Bool[], basal_extr = Bool[],
			x_A = Float64[], y_A = Float64[],
			x_B = Float64[], y_B = Float64[],
			x_P = Float64[], y_P = Float64[],
			x_R = Float64[], y_R = Float64[],
			x_init = Float64[], y_init = Float64[],
			x_emt = Float64[], y_emt = Float64[],
			x_end = Float64[], y_end = Float64[],
			t_basal_extr = Float64[])

    for i = 1:N
        type = Dict(get_type(s,i))
        label = s[i].label
        name = type["name"]
        if name != cell_type
            continue
        end

		if i > p.epi.N_init
			@error "For stochastic explorations, the 'exp' cell type is not allowed to perform cell division!"
		end

		inm_start = getvalue(type,"INM.sim_time_start")
		if isnothing(inm_start)
			inm = !maybe_no_inm
		else
			inm = !isinf(inm_start)
		end

		running_mode = s[i].running_mode
        ta = getvalue(type,"loss_apical_adhesion.sim_time_start")
        tb = getvalue(type,"loss_basal_adhesion.sim_time_start")
        tp = getvalue(type,"loss_polarity.sim_time_start")
        tr = getvalue(type,"run.sim_time_start")
        row = [name, rep, label, i, inm, ta, tb, tp, tr, running_mode, false, false]

        #row[end] = SEMTor.nuclei_is_below_basal_layer(s, p, i, basal_info)

		A_eff = SEMTor.project_onto_apical_layer(s, i)
		B_eff = SEMTor.project_onto_basal_layer(s, i)
		row[end-1] = s.X[2,i] > A_eff[2]

		xA, yA = Inf64, Inf64
		xB, yB = Inf64, Inf64
		xP, yP = Inf64, Inf64
		xR, yR = Inf64, Inf64
		xInit, yInit = Inf64, Inf64
		xEmt, yEmt = Inf64, Inf64
		xEnd, yEnd = Inf64, Inf64

		s_in = states[1]
		for i_in in 1:num_cells(states[1])
	        if cell_type == get_type(s_in, i_in).name
				tA_in = SEMTor.get_cell_event(s_in, i_in, "loss_apical_adhesion").sim_time_start
				tB_in = SEMTor.get_cell_event(s_in, i_in, "loss_basal_adhesion").sim_time_start
				tP_in = SEMTor.get_cell_event(s_in, i_in, "loss_polarity").sim_time_start
				tR_in = SEMTor.get_cell_event(s_in, i_in, "run").sim_time_start
				if (ta, tb, tp, tr) == (tA_in, tB_in, tP_in, tR_in)

					A_eff_in = SEMTor.project_onto_apical_layer(s_in, i_in)
					B_eff_in = SEMTor.project_onto_basal_layer(s_in, i_in)

					xInit = (s_in.X[1,i_in] - minimum(s_in.X[1,:])) / (maximum(s_in.X[1,:]) - minimum(s_in.X[1,:]))
					yInit = (s_in.X[2,i_in] - B_eff_in[2]) / (A_eff_in[2] - B_eff_in[2])
				end
			end
		end

		for ob in obs
			if ob.cell_sign == (ta,tb,tp)
				if ob.event == 'A'
					xA, yA = ob.rel_pos
				elseif ob.event == 'B'
					xB, yB = ob.rel_pos
				elseif ob.event == 'P'
					xP, yP = ob.rel_pos
				elseif ob.event == 'R'
					xR, yR = ob.rel_pos
				end

				if ob.t == min(ta,tb,tp,tr)
					xEmt, yEmt = ob.rel_pos
				end
			end
		end

		xEnd = (s.X[1,i] - minimum(s.X[1,:])) / (maximum(s.X[1,:]) - minimum(s.X[1,:]))
		yEnd = (s.X[2,i] - B_eff[2]) / (A_eff[2] - B_eff[2])


		x .= (state.X[2,i] for state in states)
		t_B = basal_extr_criterium(x, wnd) * p.sim.dt

		row[end] = t_B > 0

        push!(df, vcat(row,xA,yA,xB,yB,xP,yP,xR,yR,xInit,yInit,xEmt,yEmt,xEnd,yEnd,t_B))
    end
    return df
end

function sort_scenarios!(scenarios)
	sort!(scenarios, by = s -> s == "∅" ? (-1,"A") : (length(s),s))
end

function collect_times(s, p, i, event, obs)
	if (event.name in ["loss_apical_adhesion","loss_basal_adhesion", "loss_polarity", "run"])
		A_eff = SEMTor.project_onto_apical_layer(s, i)
		B_eff = SEMTor.project_onto_basal_layer(s, i)

		rel_pos_x = (s.X[1,i] - minimum(s.X[1,:])) / (maximum(s.X[1,:]) - minimum(s.X[1,:]))
		rel_pos_y = (s.X[2,i] - B_eff[2]) / (A_eff[2] - B_eff[2])

		tA = SEMTor.get_cell_event(s, i, "loss_apical_adhesion").sim_time_start
		tB = SEMTor.get_cell_event(s, i, "loss_basal_adhesion").sim_time_start
		tP = SEMTor.get_cell_event(s, i, "loss_polarity").sim_time_start
		tR = SEMTor.get_cell_event(s, i, "run").sim_time_start

		push!(obs, (event = (event.name[1] == 'r' ? 'R' : uppercase(event.name[6])),
					t = event.sim_time_start,
					rel_pos = (rel_pos_x, rel_pos_y),
					cell_sign = (tA, tB, tP, tR)))
	end
end