# a few helpers to deal with
function list_params(d::Dict, table = DataFrame(key=[],value=[]), prefix="")
	for (k,v) in d
		if typeof(v) <: Dict || typeof(v) <: Array
			list_params(v, table, prefix * k * ".")
		else
			append!(table, DataFrame(key=prefix * k, value=v))
		end
	end
	return table
end

function list_params(d::Array, table = DataFrame(key=[],value=[]), prefix="")
	for v in d
		if typeof(v) <: Dict || typeof(v) <: Array
			list_params(v, table, prefix * v["name"] * ".")
		else
			append!(table, DataFrame(key=prefix * v["name"], value=v))
		end
	end
	return table
end

function cell_events_table(d)
	cell_types = deepcopy(d["cell_types"])

	cell_params = sort(setdiff(collect(keys(cell_types[1])), ["name","cell_events","special_cell_events","life_span"]))
	append!(cell_params, ["life_span.min", "life_span.max"])

	df_p = DataFrame()
	df_ct = DataFrame(parameter = cell_params)
	df_ce = DataFrame()
	df_sce = DataFrame()

	d = deepcopy(d)
	delete!(d, "cell_types")
	p = list_params(d)
	append!( df_p, p)

	for cell_type in cell_types
		cell_type["life_span.min"] = string(cell_type["life_span"]["min"])
		cell_type["life_span.max"] = string(cell_type["life_span"]["max"])
		delete!(cell_type, "life_span")
		df_ct[!, Symbol(cell_type["name"])] .= [cell_type[p] for p in cell_params]


		cell_events = cell_type["cell_events"]

		for cell_event in cell_events
			cell_event["cell_type"] = cell_type["name"]

			df = string.(DataFrame(cell_event))
			append!(df_ce, select(df,
					[:cell_type, :name, :symbol, :factor, :abs_value, :sim_time_start, :sim_time_end, :cell_ref_time, :cell_time_start]))
		end


		spical_cell_events = cell_type["special_cell_events"]
		for cell_event in deepcopy.(spical_cell_events)
			cell_event["cell_type"] = cell_type["name"]

			df = string.(DataFrame(cell_event))

			append!(df_sce, select(df,
					[:cell_type, :name, :julia_function, :sim_time_start, :sim_time_end, :cell_ref_time, :cell_time_start]))
		end
	end

	return df_p, df_ct, df_ce, df_sce
end


function save_xlsx_parameters(fn, params::Dict; kw_args...)
	df_p, df_ct, df_ce, df_sce = cell_events_table(params)
	XLSX.writetable(fn,
		parameters=( collect(DataFrames.eachcol(df_p)), DataFrames.names(df_p) ),
		cell_type=( collect(DataFrames.eachcol(df_ct)), DataFrames.names(df_ct) ),
		cell_events=( collect(DataFrames.eachcol(df_ce)), DataFrames.names(df_ce) ),
		special_cell_events=( collect(DataFrames.eachcol(df_sce)), DataFrames.names(df_sce) ); kw_args...)
end
