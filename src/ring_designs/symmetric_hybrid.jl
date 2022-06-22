function get_symmetric_hybrid_FODO(R0=95.49, n_fodo=24)
    R0 = 95.49
    deflector_arc = 2 * pi / n_fodo / 2
    elements = RingElement[]

    G_p = 1.792847

    proton = Particle(name = "proton+",
                      m = 1.672621898e-27,
                      G = G_p,
                      gamma = sqrt(1 + 1/G_p))
    fodos = n_fodo

    quad_length = .4   # e-3
    k = 0.20488
    quad_buffer = 1.88 # (straight_section_length / (fodos*2) - quad_length)/2


    push!(elements, ElectricBendingSection(R0, proton.chiE0 / R0, deflector_arc*R0, true, proton))

    push!(elements, Drift(quad_buffer))
    push!(elements, MagneticQuadrupole(k, quad_length))
    push!(elements, Drift(quad_buffer))

    push!(elements, ElectricBendingSection(R0, proton.chiE0 / R0, deflector_arc*R0, true, proton))

    push!(elements, Drift(quad_buffer))
    push!(elements, MagneticQuadrupole(-k, quad_length))
    push!(elements, Drift(quad_buffer))

    return elements
end

function symmetric_hybrid_ring(
    ring_elements = vcat([get_symmetric_hybrid_FODO() for i in 1:24]...),
    particle = Proton() 
    )

    ring = RingStructure(ring_elements)


    p = RingParameters{Float64, eltype(ring.ringElements)}(ring=ring, particle=particle, internal_max_time_step = 1)
    return p
end
export symmetric_hybrid_ring

function get_symmetric_magnetic_FODO(R0=95.49, n_fodo=24)
    R0 = 95.49
    deflector_arc = 2 * pi / n_fodo / 2
    elements = RingElement[]

    G_p = 1.792847

    proton = Particle(name = "proton+",
                      m = 1.672621898e-27,
                      G = G_p,
                      gamma = sqrt(1 + 1/G_p))

    quad_length = .4   # e-3
    k = 0.20488
    quad_buffer = 1.88 # (straight_section_length / (fodos*2) - quad_length)/2

    rigidity = proton.momentum / e
    B0 = rigidity / R0


    push!(elements, MagneticBendingSection(R0 = R0, By = B0, length = deflector_arc*R0, particle = proton, n=0.0))

    push!(elements, Drift(quad_buffer))
    push!(elements, MagneticQuadrupole(k, quad_length))
    push!(elements, Drift(quad_buffer))

    push!(elements, MagneticBendingSection(R0 = R0, By = B0, length = deflector_arc*R0, particle = proton, n=0.0))

    push!(elements, Drift(quad_buffer))
    push!(elements, MagneticQuadrupole(-k, quad_length))
    push!(elements, Drift(quad_buffer))

    return elements
end

function symmetric_magnetic_ring(
    ring_elements = vcat([get_symmetric_magnetic_FODO() for i in 1:24]...),
    particle = Proton() 
    )

    ring = RingStructure(ring_elements)
    p = RingParameters{Float64, eltype(ring.ringElements)}(ring=ring, particle=particle, internal_max_time_step = 1)
    return p
end
export symmetric_magnetic_ring
