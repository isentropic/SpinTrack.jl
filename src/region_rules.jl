function change_region!(integrator, cw)
    this_element_index = integrator.p.element_index

    N = integrator.p.ring.length_of_ring_elements
    if cw
        next_ring_element_index = prev_element_index(this_element_index, N)
    else
        next_ring_element_index = next_element_index(this_element_index, N)
    end
    next_ring_element = integrator.p.ring.ringElements[next_ring_element_index]
    integrator.opts.dtmax = next_ring_element.length
    integrator.p.element_index = next_ring_element_index
    return nothing
end

function next_region_time(integrator)
    ringElement = integrator.p.ring.ringElements[integrator.p.element_index]
    return ringElement.length + integrator.t
end

function increaseLongitudinalEnergy!(integrator, ΔW)
    integrator.u[6] += ΔW/integrator.p.particle.K0
end


function next_rf_time(integrator, cw)
    ring = integrator.p.ring
    return integrator.t  + ring.total_length
end

function apply_rf!(integrator, cw)
    time = integrator.u[time_index]
    V = integrator.p.RF_voltage

    ΔW = e * V * sin(integrator.p.RF_omega * time + integrator.p.RF_phase)
    increaseLongitudinalEnergy!(integrator, ΔW)
end
