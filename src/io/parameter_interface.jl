function load_parameters(fn)
	param = nothing
	sim_name_fn = ""
	if isfile(fn)
		if occursin('$', fn)
			return nothing
		end
		if lowercase(fn[end-4:end]) == ".toml"
			sim_name_fn = basename(fn)[1:end-5]
			file_str = read(fn,String)
			if file_str[1] == '\ufeff'
				file_str = file_str[4:end]
			end
			param = TOML.parse(file_str)
		end

		if lowercase(fn[end-4:end]) == ".xlsx"
			sim_name_fn = basename(fn)[1:end-5]
			param = load_xlsx_parameters(fn)
		end
	else
		@error "Cannot load parameterfile $sim_name_fn. File does not exist."
	end

	if ismissing(param)
		param = nothing
	end

	if !isnothing(param)
		sim_name = getvalue(param,"sim.name")
		if ismissing(sim_name) || isempty(sim_name)
			setvalue!(param,"sim.name",sim_name_fn)
		end
	end

	return param
end

function setvalue!(d::Dict, path::String, val, insert=false)
    idx = findfirst( '.', path)
    if isnothing(idx)
        d[path] = val
    else
        base = path[1:idx-1]
        rest = path[idx+1:end]

        if haskey(d, base)
            setvalue!(d[base], rest, val, insert)
        elseif insert
            d[base] = Dict{String,Any}()
            setvalue!(d[base], rest, val, insert)
        else
            # it could be either a cell_event/special_cell_event or a cell_type
            for array_key = ["cell_types","cell_events","special_cell_events"]
                if haskey(d, array_key)
                    setvalue!(d[array_key], path, val)
                end
            end
        end
    end
end

setvalue!(d, ::Missing, ::Missing, ::Bool) = nothing

function setvalue!(d::Vector, path::String, val)
    idx = findfirst( '.', path)
    if isnothing(idx)
        @warn "unexpected parameter structure"
    else
        base = path[1:idx-1]
        rest = path[idx+1:end]
        for x in d
            if x["name"] == base
                setvalue!(x, rest, val)
            end
        end
    end
end


function getvalue(d::Dict, path::String)

	if haskey(d, path)
		return d[path]
	end

    idx = findfirst( '.', path)
	if isnothing(idx)
		idx = 0
	end

    base = path[1:idx-1]
    rest = path[idx+1:end]

    if haskey(d, base)
        return getvalue(d[base], rest)
    else
        # it could be either a cell_event/special_cell_event or a cell_type
        for array_key = ["cell_types","cell_events","special_cell_events"]
            if haskey(d, array_key)
                v = getvalue(d[array_key], path)
				if !isnothing(v)
					return v
				end
            end
        end
    end

	return nothing
end
function getvalue(d::Vector, path::String)

    idx = findfirst( '.', path)
	if isnothing(idx)
		idx = length(path)+1
	end


    base = path[1:idx-1]
    rest = path[idx+1:end]
    for x in d
        if x["name"] == path
            return x
        end
        if x["name"] == base
            return getvalue(x, rest)
        end
    end

	return nothing
end



function flatten(d::Dict; prefix = "", p = Dict())
    for (key, value) in d
        if isa(value, Dict)
            flatten(value; prefix=prefix * key * ".", p)
        elseif isa(value, Array)
            flatten(value; prefix=prefix, p)
        else
            p[prefix * key] = value
        end
    end
    return p
end


function flatten(a::Array; prefix = "", p = Dict())
    for value in a
        if isa(value, Dict) && haskey(value, "name")
            flatten(value; prefix=prefix * value["name"] * ".", p)
        end
    end
    return p
end

function Dict(ct::EpiCellType)
	p = Dict{String,Any}(string(fn) => getfield(ct,fn) for fn in fieldnames(EpiCellType) if fn ∉ [:cell_events, :special_cell_events, :life_span])

	p["life_span.min"] = ct.life_span.min
	p["life_span.max"] = ct.life_span.max

	for ce in ct.cell_events
		ce_name = ce.name
		for fn in fieldnames(CellEvent)
			p[ce_name * "." * string(fn)] = getfield(ce, fn)
		end
	end

	for ce in ct.special_cell_events
		ce_name = ce.name
		for fn in fieldnames(SpecialCellEvent)
			if fn != :value
				p[ce_name * "." * string(fn)] = getfield(ce, fn)
			end
		end
	end
	return p
end
