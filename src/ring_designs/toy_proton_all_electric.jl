function toy_all_electric(particle, n; sections=3, kwargs...)
    R0 = 100
    rigidity = particle.momentum / e
    Ex = rigidity * particle.beta * c / R0



    ringElements = fill(ElectricBendingSection(R0, Ex, 2pi * R0/sections, true, particle, n + 1, false, 0.0), sections)
    ring = RingStructure(ringElements)
    return RingParameters{Float64, eltype(ring.ringElements)}(particle=particle, ring=ring; kwargs...)
end

function toy_proton_all_electric(n=0.1; sections=3, kwargs...)
    return toy_all_electric(Proton(), n; sections=sections, kwargs...)
end
export toy_proton_all_electric

function toy_muon_all_electric(n=0.1; sections=3, kwargs...)
    return toy_all_electric(Muon(), n; sections=sections, kwargs...)
end
export toy_muon_all_electric


function toy_all_magnetic(particle, n; sections=3, kwargs...)
    R0 = 100
    rigidity = particle.momentum / e
    B0 = rigidity / R0



    ringElements = fill(MagneticBendingSection(R0 = R0, By = B0, length = 2pi * R0/sections, particle = particle, n=n),
                        sections)
    ring = RingStructure(ringElements)
    return RingParameters{Float64, eltype(ring.ringElements)}(particle=particle, ring=ring; kwargs...)

end

function toy_proton_all_magnetic(n=0.1; sections=3, kwargs...)
    return toy_all_magnetic(Proton(), n, sections=sections, kwargs...)
end
export toy_proton_all_magnetic

function toy_muon_all_magnetic(n=0.1; sections=3, kwargs...)
    return toy_all_magnetic(Muon(), n; sections=sections, kwargs...)
end
export toy_muon_all_magnetic
