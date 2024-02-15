module SEMTor

using UnPack
using Printf
using LinearAlgebra
using Distributions
using Random
using StructArrays
using StaticArrays, HybridArrays

# config files
using DataFrames     # takes time
using TOML
using XLSX           # takes time

using GLMakie
using CairoMakie

import Base.show
include("core/definitions.jl")
include("cell_events/cell_events.jl")

include("core/edge_list.jl")
include("core/verlet_list.jl")

include("core/cell_type.jl")
include("io/parameter_interface.jl")
include("core/random_params.jl")
include("core/parameters.jl")
include("core/state.jl")
include("cell_events/cell_event_functions.jl")
include("core/init_state.jl")

include("algorithms/pbd.jl")
include("core/compute_forces.jl")
include("core/simulate.jl")

include("cell_events/cell_division.jl")
include("cell_events/loss_apical_basal_connection.jl")
include("cell_events/reset_cell_cycle.jl")

include("visualisation/plot_makie.jl")
include("analysis/compute_statistics.jl")

include("io/xlsx_export.jl")
include("io/load_xlsx_parameters.jl")
include("io/save_xlsx_parameters.jl")

export CellPhase, CellInitTypes, ReferenceTimePoint
export PBD
export TimeInterval, Rectangle, AffineTraffo
export State, State2D, State3D, StaticParameters
export get_type, by_label

export divide_cell, loss_apical_connection, loss_basal_connection
export simulate, simulate_ensemble

export read_json, write_json, load_json, save_json

export EpiCell, ApicalCytoskeleton, BasalCytoskeleton, ApicalJunction, BasalJunction, BasalCytoskeleton, ApicalParticle, BasalParticle
export CellEvent, SpecialCellEvent
export divide_cell, isCellCycleEvent, loss_apical_connection, loss_basal_connection
export create_cell_dataframe, create_connection_dataframe, df_to_xlsx
export load_state, save_state
export bio_connections
export minimal_plot_rect
export plot_state!, plot_state
export export_statistics_to_xlsx

export get_age, get_cell_phase, G1, G2, mitosis, num_cells, get_cell_phase_idx

export load_xlsx_parameters, save_xlsx_parameters, load_parameters

export compute_statistics

export getvalue, setvalue!
export quickcopy
export load_xlsx_parameters, save_xlsx_parameters


# include("precompile_includer.jl")
end
