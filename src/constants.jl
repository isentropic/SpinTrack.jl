
getFields(u, p, t, ringElement::RingElement) = (zeros(SVector{3}), zeros(SVector{3}), 0.0)
getCurvature(ringElement)::Float64 = 0.0

struct RingStructure{ElementType}
    ringElements::Vector{ElementType}
    total_length::Float64
    length_of_ring_elements::Int
    function RingStructure(ringElements)
        united = Union{unique(typeof.(ringElements))...}
        return new{united}(ringElements, sum(e.length for e in ringElements), length(ringElements))
    end
end

@with_kw mutable struct RingParameters{R<:Real,ElementType<:RingElement}
    @deftype R

    particle::Particle = Particle()
    ring::RingStructure{ElementType} = getRingStructure()
    region_change_function!::Function = change_region!

    # Simulation related values (dynamic)
    element_index::Int = 1
    verbose::Bool = false

    alg::Any = Tsit5()

    save_positions::Bool = false

    starting_time = 0.0
    turns = 500.0

    saving_time_step = ring.total_length / 10
    internal_max_time_step = ring.ringElements[element_index].length

    state_tolerance_rel::Union{Float64,Vector{Float64}} = 1e-12
    state_tolerance_abs::Union{Float64,Vector{Float64}} = 1e-14 # reduce abs tolerance to get long term running

    # EDM
    #
    if_EDM_on::Bool = false
    η = 1.9e-15

    # RF
    RF_on::Bool = false
    RF_phase = 0.0
    cyclotron_freq = 1 / (ring.total_length / particle.beta / c)
    RF_omega = 2pi * cyclotron_freq * 80   # 30 times faster than
    RF_voltage = 1.89e3

    is_losing_energy::Bool = false
    beam_current = particle.beta * c * e / ring.total_length
    ring_impedance = 1e3
    energy_loss_per_meter = beam_current^2 * ring_impedance / ring.total_length

    global_E_y::Bool = false
    E_y = 0.0

    global_B_R::Bool = false
    B_multipole = 0
    B_R = 1e-14  # 1pT B field
    B_R_phase = 0.0

    global_B_L::Bool = false
    B_L = 1e-9  # 1nT B field

    global_B_y::Bool = false
    B_y = 1e-12
end

next_element_index(i, N) = (i) % (N) + 1
prev_element_index(i, N) = (i + (N - 2)) % (N) + 1
