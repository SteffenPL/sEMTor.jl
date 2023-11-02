function reset_cell_cycle(s::State{d}, p::StaticParameters, i::Int64) where {d}

    cell1 = new_cell_copy_startimes(s, p, i)
    cell1.age = 0.0
    cell1.birth_time = s.t
    type = get_type(s,i)

    if type.life_span.min < type.life_span.max
        cell1.division_time = rand(Uniform(type.life_span.min, type.life_span.max))
    elseif type.life_span.min ≈ type.life_span.max
        cell1.division_time = type.life_span.min
    elseif !isinf(type.life_span.min)
        @error( "life_span.min > life_span.max" )
    end

    s[i] = cell1
    type = get_type(s,i)

    for event in type.cell_events
        init_cell_event!(s, p, i, event)
    end
    for event in type.special_cell_events
        init_cell_event!(s, p, i, event)
    end
end
