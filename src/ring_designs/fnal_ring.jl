function FNAL_ring(n=0.1)
    R0 = 7.112
    d = Deuteron()

    rigidity = d.momentum / e
    By = rigidity / R0 / (1 + d.G  * c * d.beta * d.gamma^2 / (1 - d.G * d.beta^2 * d.gamma^2)/(d.beta * c))


    Ex =  d.G * By * c * d.beta * d.gamma^2 / (1 - d.G * d.beta^2 * d.gamma^2)
#     @show By, Ex

    ringElements = [HybridBendingSection(R0, By, Ex, n, 2pi*R0)]
    ring = RingStructure(ringElements)
end

FNAL_deutron_ring(n=0.1) = RingParameters{Float64, eltype(FNAL_ring(n).ringElements)}(particle=Deuteron(), ring=FNAL_ring(n));
export FNAL_deutron_ring

function FNAL_muon_ring(k = 17e6 * .43; kwargs...)
    R0 = 7.112
    d = Muon()
    rigidity = d.momentum / e
    B0 = rigidity / R0
    total_length = 2pi * R0
    quad_length = total_length * 0.43 / 4
    drift_length = total_length / 4  - quad_length
    # 43 %
    fodo = [CurvedElectricQuadrupole(R0 = R0,
                                     length = quad_length,
                                     k = k;
                                     kwargs...),
            CurvedDrift(R0, drift_length)]
    ring = RingStructure(repeat(fodo, 4))

    return RingParameters{Float64, eltype(ring.ringElements)}(particle=d, ring=ring, global_B_y = true, B_y = B0; )
end
export FNAL_muon_ring
# Muon
# const g = 2.0023318418;
# const a = (g-2)/2
# const m = 1.883531475e-28; //# mass of a muon (kgs)
