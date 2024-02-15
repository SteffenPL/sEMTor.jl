function apply!(ce::CellEvent, s::State, p::StaticParameters, i::Int64)
    setproperty!(s[i], ce.symbol,
                ce.abs_value
                + ce.factor * getproperty(s[i], ce.symbol))
end

function apply!(ce::SpecialCellEvent, s::State, p::StaticParameters, i::Int64)
    ce.julia_function(s, p, i)
end

function get_relative_time(s::State, i::Int64, ce::AbstractCellEvent)
    ce_age = get_age(s, i)

    if ce.cell_ref_time == Cell_G2_Start
        ce_age -= (s.cells.division_time[i]
                - s.cells.duration_mitosis[i]
                - s.cells.duration_G2[i])
    elseif ce.cell_ref_time == Cell_Mitosis_Start
        ce_age -= (s.cells.division_time[i]
                - s.cells.duration_mitosis[i])
    elseif ce.cell_ref_time == Cell_Division
        ce_age -= s.cells.division_time[i]
    end

    return ce_age
end

function is_active(ce, t, cell_age)
    return (ce.sim_time_start <= t < ce.sim_time_end) && (ce.cell_time_start <= cell_age)
end


function is_starting(ce, t::Float64, cell_age::Float64, dt::Float64)
    return (ce.sim_time_start <= t < ce.sim_time_start+dt) && (ce.cell_time_start <= cell_age) ||
           (ce.sim_time_start < t < ce.sim_time_end) && (ce.cell_time_start <= cell_age < ce.cell_time_start + dt)
end

function update!(s::State, p::StaticParameters, i::Int64, ce::AbstractCellEvent)
    rel_time = get_relative_time(s, i, ce)

    if is_starting(ce, s.t, rel_time, p.alg.dt)
        apply!(ce, s, p, i)
        return true
    end
    return false
end

function init_cell_event!(state::State, p::StaticParameters, i::Int64, ce::AbstractCellEvent)
    rel_time = get_relative_time(state, i, ce)

    if is_active(ce, state.t, rel_time)
        apply!(ce, state, p, i)
    end
end
