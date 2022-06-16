const c = 2.99792458e8
const e = 1.6021766208e-19
export c, e

@with_kw struct Particle @deftype Float64
    name::String = "proton+"
    m = 1.672621898e-27
    G = 1.792847
    gamma = sqrt(1 + 1/G)
    beta = sqrt(1 - (1/gamma)^2)
    momentum = m * beta * c * gamma
    Energy0 = sqrt((momentum*c)^2 + (m*c^2)^2)
    eta0 = gamma - 1.0
    K0 = eta0 * m * c^2

    chiM0 = momentum / e
    chiE0 = momentum * beta * c / e
    charge = 1.0 # measured in e
    kappa = -beta * c * gamma / (gamma+1)
end



Proton() = Particle(name = "proton+",
                    m = 1.672621898e-27,
                    G = 1.792847,
                    gamma = sqrt(1 + 1/1.792847))

Muon() = Particle(name = "muon+",
                  m = 1.883531627e-28,
                  G = 0.00116591810)

Deuteron() = Particle(name = "deuteron+",
                    m = 1.672621898e-27 * 2,
                    G = -0.1425617692,
                    gamma = 1/sqrt(1 - 0.2574646862654724^2))
