function read_cell_event(table, event_name="cell event")

    events = Dict{String,Any}()
    for row in 2:size(table,1)
        event = Dict{String,Any}()
        ct_name = table[row,1]

        if ismissing(ct_name)
            continue
        end

        for col in 2:size(table,2)
            key = table[1,col]
            val = table[row,col]
            if !ismissing(key)
                if !ismissing(val)
                    setvalue!(event, key, val, true)
                else
                    @warn "A '$(event_name)' for cell type '$(ct_name)' has a missing value for key '$(key)'."
                end
            end
        end
        if !haskey(events, ct_name)
            events[ct_name] = []
        end
        push!(events[ct_name], event)
    end
    return events
end

function load_xlsx_parameters(fn; hide_warning=false)
    XLSX.openxlsx(fn, mode="r") do xf
        if !("parameters" in XLSX.sheetnames(xf))
            if !hide_warning
                @warn "The parameter file $fn does not contain a sheet called 'parameters'."
            end
            return nothing
        end

        parameters = xf["parameters"][:]
        cell_types = xf["cell_type"][:]
        cell_events = xf["cell_events"][:]
        special_cell_events = xf["special_cell_events"][:]

        p = Dict{String,Any}()
        for row in 2:size(parameters,1)
            key = parameters[row,1]
            val = parameters[row,2]
            setvalue!(p, key, val, true)
        end

        events = read_cell_event(cell_events,"cell event")
        special_events = read_cell_event(special_cell_events,"special cell event")

        p["cell_types"] = []
        for col in 2:size(cell_types,2)
            ct_name = cell_types[1,col]
            ct_dict = Dict{String,Any}()
            ct_dict["name"] = ct_name

            if ismissing(ct_name)
                continue
            end

            for row in 2:size(cell_types,1)
                key = cell_types[row,1]
                value = cell_types[row,col]
                if ismissing(key)
                    continue
                end
                if ismissing(value)
                    @warn "Not value for parameter '$(key)' of cell type '$(ct_name)'."
                    continue
                end
                setvalue!(ct_dict, key, value, true)
            end
            ct_dict["cell_events"] = (haskey(events, ct_name) ? events[ct_name] : [])
            ct_dict["special_cell_events"] =(haskey(events, ct_name) ? special_events[ct_name] : [])
            push!(p["cell_types"], ct_dict)
        end
        return p
    end
end
