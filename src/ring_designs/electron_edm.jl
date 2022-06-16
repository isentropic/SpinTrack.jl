function change_region_electron_edm!(integrator, cw)
    ringElement = integrator.p.ring.ringElements[integrator.u.ringElementIndex]

    if cw
        nextRingElementIndex = getPrevRingElementIndex(integrator.u.ringElementIndex, integrator.p.ring.length_of_ring_elements)
    else
        nextRingElementIndex = getNextRingElementIndex(integrator.u.ringElementIndex, integrator.p.ring.length_of_ring_elements)
    end

    nextRingElement = integrator.p.ring.ringElements[nextRingElementIndex]
    orientation = ringElement.curvature
    next_orientation = nextRingElement.curvature
    integrator.p.particle = nextRingElement.particle

    if orientation != next_orientation
        integrator.u[1] *= -1
        integrator.u[2] *= -1
        integrator.u[7] *= -1
    end
    for c in full_cache(integrator)
        c.ringElementIndex = nextRingElementIndex
        c.elementEnteringS = integrator.t
    end
end


function electron_EDM_ring()
    e1 = Particle(
        name = "electron-",
        m = 9.10938356e-31,
        G = 0.001159652181643,
       # gamma = sqrt(1 + 1/0.001159652181643),
        gamma=1.4
    )

    e2 = Particle(
        name = "electron-",
        m = 9.10938356e-31,
        G = 0.001159652181643,
        #gamma = 15.0,
        gamma=2.6
    )
    R1 = 0.07
    R2 = 0.23

    ring_elements = [
        ElectricBendingSection(R1, e1.chiE0 / R1, pi * R1, 1, e1),
        ElectricBendingSection(R2, e2.chiE0 / R2, pi * R2, 1, e2),
        ElectricBendingSection(R1, e1.chiE0 / R1, pi * R1, 0, e1),
        ElectricBendingSection(R2, e2.chiE0 / R2, pi * R2, 0, e2)
    ]
    ring = RingStructure(ring_elements)

    p = RingParameters{Float64, eltype(ring.ringElements)}(ring=ring, particle=e1,);
    p.saving_time_step = 0.32
    p.vacuum_max_step_size = 5e-10
    p.region_change_function! = change_region_electron_edm!
    return p
end
export electron_EDM_ring
