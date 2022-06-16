module SpinTrack

using DiffEqCallbacks
using DSP
using FFTW
using Interpolations
using LinearAlgebra
using LsqFit
using Measurements
using OrdinaryDiffEq
using Parameters
using RecipesBase
using Reexport
using StaticArrays
using Statistics

@reexport using LaTeXStrings

const M = Measurements
const time_index = 5
const size_u = 9

u0_long = zeros(size_u); u0_long[9] = 1.0
u0_rad = zeros(size_u); u0_rad[7] = 1.0
u0_vert = zeros(size_u); u0_vert[8] = 1.0

u1_long(oscillation_amplitude) = [oscillation_amplitude, 0, oscillation_amplitude, 0, 0, 0, 0, 0, 1, ]
u1_rad(oscillation_amplitude)  = [oscillation_amplitude, 0, oscillation_amplitude, 0, 0, 0, 1, 0, 0, ]
u1_vert(oscillation_amplitude) = [oscillation_amplitude, 0, oscillation_amplitude, 0, 0, 0, 0, 1, 0, ]

u_x(oscillation_amplitude) = [oscillation_amplitude, 0, 0, 0, 0, 0, 0, 0, 1, ]
u_y(oscillation_amplitude) = [0, 0, oscillation_amplitude, 0, 0, 0, 0, 0, 1, ]
u_p(offset, particle) = [0, 0, 0.0, 0, 0, get_delta(offset, particle), 0, 0, 1, ]

export u_x, u_y, u_p
export u0_long, u0_rad, u0_vert, u1_long, u1_rad, u1_vert

include("particles.jl")
include("ring_designs/base_elements.jl")
include("constants.jl")
include("ring_designs/fnal_ring.jl")
include("ring_designs/electron_edm.jl");
include("ring_designs/symmetric_hybrid.jl");
include("ring_designs/toy_proton_all_electric.jl");

include("global_em_fields.jl")
include("region_rules.jl")
include("bmt_equation.jl")

include("solutions.jl")
include("analysis/analysis.jl")
include("analysis/plottings.jl")
include("utils.jl")

export get_solution, get_moving_average, get_simple_linear_rate

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end
end
