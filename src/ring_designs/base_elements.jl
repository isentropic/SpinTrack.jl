abstract type RingElement end
struct Drift <: RingElement
    length::Float64
end
function getFields(u, p, t, element::Drift)::Tuple{StaticArrays.SVector{3, Float64}, StaticArrays.SVector{3, Float64}, Float64}
    zeros(SVector{3}), zeros(SVector{3}), 0.0
end

function get_fields!(ebv, u, p, t, element::Drift)
    return nothing
end

struct ElectricBendingSection <: RingElement
    R0::Float64
    Ex::Float64
    length::Float64
    curvature::Bool
    particle::Particle

    n::Float64

    is_Ey_compensated::Bool
    Ey::Float64
end

getCurvature(ringElement::ElectricBendingSection)::Float64 = 1/ringElement.R0

ElectricBendingSection(R0, Ex, length, curvature, particle) = ElectricBendingSection(R0, Ex, length, curvature, particle, 1.0, false, 0.0)

function getFields(u, p, t, element::ElectricBendingSection)::Tuple{StaticArrays.SVector{3, Float64}, StaticArrays.SVector{3, Float64}, Float64}
    x, y = u[1], u[3]

    E0 = element.Ex
    R0 = element.R0

    n = element.n
    E_R = - E0 * (1 - n*x/R0 + n*(n+1) * x^2/R0^2/2)
    V = E0 * (x - n*x^2/(2*R0) - x^3 * (n*(n+1)/6/R0^2))  # Taylor expansion of the log

    E_z = - E0 * ((n-1) * y / (R0))
    V += E0 * ((n-1)*y^2/(2*R0))

    if element.is_Ey_compensated
        extra_Ez = element.Ey
        E_z += extra_Ez

        V += - extra_Ez * y
    end

    Efield = SA[E_R, E_z, 0.0]
    Bfield = SA[0.0, 0.0, 0.0]
    return Efield, Bfield, V
end

function get_fields!(ebv, u, p, t, element::ElectricBendingSection)
    x, y = u[1], u[3]

    E0 = element.Ex
    R0 = element.R0

    n = element.n
    E_R = - E0 * (1 - n*x/R0 + n*(n+1) * x^2/R0^2/2)
    V = E0 * (x - n*x^2/(2*R0) - x^3 * (n*(n+1)/6/R0^2))  # Taylor expansion of the log

    E_z = - E0 * ((n-1) * y / (R0))
    V += E0 * ((n-1)*y^2/(2*R0))

    if element.is_Ey_compensated
        extra_Ez = element.Ey
        E_z += extra_Ez

        V += - extra_Ez * y
    end

    ebv[1] += E_R
    ebv[2] += E_z
    ebv[7] += V

    return nothing
end

@with_kw struct MagneticBendingSection <: RingElement @deftype Float64
    R0
    By
    length
    particle::Particle
    n
end

getCurvature(ringElement::MagneticBendingSection)::Float64 = 1/ringElement.R0

function getFields(u, p, t, element::MagneticBendingSection)::Tuple{StaticArrays.SVector{3, Float64}, StaticArrays.SVector{3, Float64}, Float64}
    x, y = u[1], u[3]

    B0 = element.By
    R0 = element.R0
    n = element.n
    By = B0 - n * B0 /R0 * x # + n * B0 / R0 * y^2 / 2 / R0
    Bx = - n * B0 /R0 * y
    Efield = SA[0.0, 0.0, 0.0]
    Bfield = SA[Bx, By, 0.0]
    V = 0.0
    return Efield, Bfield, V
end


struct HybridBendingSection <: RingElement
    R0::Float64
    By::Float64
    Ex::Float64

    n::Float64
    length::Float64
end

getCurvature(ringElement::HybridBendingSection)::Float64 = 1/ringElement.R0
function getFields(u, p, t, element::HybridBendingSection)::Tuple{StaticArrays.SVector{3, Float64}, StaticArrays.SVector{3, Float64}, Float64}
    x, y = u[1], u[3]

    E0 = -element.Ex
    R0 = element.R0

    E_R = E0 * (1-x/R0 + x^2/R0^2) # 2nd order is critical
    V = -E0 * (x - x^2/(2*R0) + x^3/3/R0^2)
    Efield = SA[E_R, 0.0, 0.0]

    B0 = element.By
    n = element.n

    Bx = -n * B0 / R0 * y
    By = B0 - n * B0/R0 * x # @TODO second order term
    Bfield = SA[Bx, By, 0.0]

    return Efield, Bfield, V
end


struct MagneticQuadrupole <: RingElement
    k::Float64
    Δx::Float64
    Δy::Float64

    length::Float64
end
MagneticQuadrupole(k, length) = MagneticQuadrupole(k, 0, 0, length)

function getFields(u, p, t, element::MagneticQuadrupole)::Tuple{StaticArrays.SVector{3, Float64}, StaticArrays.SVector{3, Float64}, Float64}
    x = u[1]
    y = u[3]

    Bx = element.k * y
    By = element.k * x

    Bfield = SA[Bx, By, 0]

    return zeros(SVector{3}), Bfield, 0.0
end
function get_fields!(ebv, u, p, t, element::MagneticQuadrupole)
    x = u[1]
    y = u[3]

    Bx = element.k * y
    By = element.k * x

    ebv[4] += Bx
    ebv[5] += By

    return nothing
end
struct CurvedDrift <: RingElement
    R0::Float64
    length::Float64
end
getCurvature(ringElement::CurvedDrift)::Float64 = 1/ringElement.R0

@with_kw struct CurvedElectricQuadrupole <: RingElement @deftype Float64
    R0
    length
    k

    isRF::Bool = false

    RF_omega = 0.0
    RF_phase = 0.0
    RF_Ey = 0.0
    RF_Br = 0.0
end

getCurvature(ringElement::CurvedElectricQuadrupole)::Float64 = 1/ringElement.R0

function getFields(u, p, t, element::CurvedElectricQuadrupole)::Tuple{StaticArrays.SVector{3, Float64}, StaticArrays.SVector{3, Float64}, Float64}
    x, y = u[1], u[3]

    k = element.k

    Ex = + k * x
    Ey = - k * y
    V = - k * x^2 / 2 + k * y^2 /2
    Efield = SA[Ex, Ey, 0]
    Bfield = zeros(SVector{3})

    if element.isRF
        Ey = element.RF_Ey
        Br = element.RF_Br
        ω = element.RF_omega
        ϕ = element.RF_phase
        t = u[5]
        Erf = SA[0, Ey * sin(ω * t + ϕ), 0]
        Brf = SA[Br * sin(ω * t + ϕ), 0, 0]
        Vrf =  - Ey * sin(ω * t + ϕ) * y

        return Efield + Erf, Bfield + Brf, V + Vrf
    end

    return Efield, Bfield, V
end
