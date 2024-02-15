Base.@kwdef mutable struct EpiCellType
    name::String = "control"
    color::String = "none"
    prototype_idx::Int64 = 1

    cell_events::Vector{CellEvent} = [CellEvent(:R_hard, 0.7, 0.0; cell_ref_time = Cell_G2_Start)]
    special_cell_events::Vector{SpecialCellEvent} = [SpecialCellEvent(divide_cell; cell_ref_time = Cell_Division)]

    # initial radius parameters
    R_soft::Float64 = 1.
    R_hard::Float64 = 0.3

    # stiffness parameters
    stiffness_repulsion::Float64 = 1.
    stiffness_nuclei_apical::Float64 = 2.
    stiffness_nuclei_basal::Float64 = 2.
    stiffness_apical_apical::Float64 = 5.
    stiffness_straightness::Float64 = 15.

    cytoskeleton_init::Float64 = 1.5
    k_cytoskeleton::Float64 = 5.

    apical_cytos_strain::Float64 = 0.
    basal_cytos_strain::Float64 = 0.
    running_mode::Int64 = 0
    running_speed::Float64 = 0.

    apical_junction_init::Float64 = 10. / 30.
    k_apical_junction::Float64 = 1.
    apical_junction_strain::Float64 = -1.

    diffusion::Float64 = 0.1  # diffusion coefficient
    basal_damping_ratio::Float64 = 1.
    max_basal_junction_dist::Float64 = 1. / 3.

    life_span::TimeInterval = TimeInterval(10.,21.)
    duration_G2::Float64 = 0.5
    duration_mitosis::Float64 = 0.5
end


function Base.show(io::IO, t::EpiCellType)
    compact = get(io, :compact, false)
    if !compact
        println(io, "Cell type ", t.name, " #", t.prototype_idx)
        for event in t.cell_events 
            println(io, event)
        end
        for event in t.special_cell_events 
            println(io, event)
        end
    else
        print(io, "Cell type ", t.name)
    end
end




Base.@kwdef mutable struct EpiCell
    type_idx::Int64
    prototype_idx::Int64

    R_soft::Float64 = 1.
    R_hard::Float64 = 0.3

    # stiffness parameters
    stiffness_repulsion::Float64 = 1.
    stiffness_nuclei_apical::Float64 = 2.
    stiffness_nuclei_basal::Float64 = 2.
    stiffness_apical_apical::Float64 = 5.
    stiffness_straightness::Float64 = 15.

    k_cytoskeleton::Float64 = 5.

    apical_cytos_strain::Float64 = 0.
    basal_cytos_strain::Float64 = 0.
    running_mode::Int64 = 0
    running_speed::Float64 = 0.

    k_apical_junction::Float64 = 1.
    apical_junction_strain::Float64 = -1.

    diffusion::Float64 = 0.1  # diffusion coefficient

    basal_damping_ratio::Float64 = NaN64
    max_basal_junction_dist::Float64 = NaN64

    duration_G2::Float64 = 0.5
    duration_mitosis::Float64 = 0.5

    birth_time::Float64 = -Inf64
    division_time::Float64 = Inf64

    has_apical_cytos::Bool = true
    has_basal_cytos::Bool = true

    apical_cytos_rest_length::Float64 = 0.
    basal_cytos_rest_length::Float64 = 0.

    age::Float64 = 0.
    #phase::CellPhase = G0
    is_running::Bool = false
    is_below::Bool = false
    on_basal_layer::Bool = true

    label::Int64 # unique label of the cell
end

const type_cell_properties = intersect(fieldnames(EpiCell), fieldnames(EpiCellType))
get_age(cell) = cell.age

function EpiCell(type::EpiCellType; kw...)
    EpiCell(;kw..., zip(type_cell_properties,
        map(prob -> getproperty(type, prob), type_cell_properties))...)
end

function init_cell!(cell, type::EpiCellType; t=-1.0)

    cell.apical_cytos_rest_length = type.cytoskeleton_init
    cell.basal_cytos_rest_length = type.cytoskeleton_init
    if isnan( cell.max_basal_junction_dist )
        cell.max_basal_junction_dist = p.epi.max_basal_junction_dist
    end
    if isnan( cell.basal_damping_ratio )
        cell.basal_damping_ratio = p.epi.basal_damping_ratio
    end

    if type.life_span.min < type.life_span.max
        cell.division_time = rand(Uniform(type.life_span.min, type.life_span.max))
    elseif type.life_span.min ≈ type.life_span.max
        cell.division_time = type.life_span.min
    elseif !isinf(type.life_span.min)
        @error( "life_span.min > life_span.max" )
    end

    if t < 0.0
        cell.birth_time = -rand(Uniform([0., cell.division_time]...))
        cell.age = -cell.birth_time
    else
        cell.birth_time = t
        cell.age = 0.0
    end

    return cell
end
