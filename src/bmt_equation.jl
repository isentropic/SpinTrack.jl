function bmt_equation(du, u, p, t; cw=true)
    m_p = p.particle.m
    q = e * p.particle.charge
    G_p = p.particle.G
    eta0 = p.particle.eta0
    chiM0 = p.particle.chiM0
    chiE0 = p.particle.chiE0
    kappa = p.particle.kappa
    p_magic = p.particle.momentum

    @fastmath begin
        ringElement = p.ring.ringElements[p.element_index]  # Type stability is crucial
        h = getCurvature(ringElement)
        if cw
            Efield, Bfield, V = get_total_fields(u, p, t, ringElement, cw=true)

            Ex, Ey, Es = Efield[1], Efield[2], Efield[3]
            Bx, By, Bs = Bfield[1], Bfield[2], Bfield[3]
        else
            u[3] = -u[3]
            Efield, Bfield, V = get_total_fields(u, p, t, ringElement, cw=false)
            u[3] = -u[3]

            Ex, Ey, Es = Efield[1], -Efield[2], -Efield[3]
            Bx, By, Bs = Bfield[1], -Bfield[2], -Bfield[3]

            Efield = SA[Ex, Ey, Es]
            Bfield = SA[Bx, By, Bs]
        end

        # Rate of position
        # Berz Introduction to Beam Physics
        # u[1] x
        # u[2] a = px/p0
        # u[3] y
        # u[4] b = py/p0
        # u[5] l = k(t-t0)
        # u[6] δ = (K-K0)/K0
        # h = h(s) = dΘ(s)/ds  CURVATURE!

        let (x, a, b, δ) = (u[1], u[2], u[4], u[6])
            eta = eta0 * (1 + δ) - e * V / m_p / c^2
            ps_over_p0 = sqrt(eta * (2 + eta) / eta0 / (2 + eta0) - a^2 - b^2)

            du[1] = a * (1 + h * x) / ps_over_p0

            first = (1 + eta) / (1 + eta0) * Ex / chiE0 / ps_over_p0 + b * Bs / chiM0 / ps_over_p0 - By / chiM0
            du[2] = (1 + h * x) * (first) + h * ps_over_p0

            du[3] = b * (1 + h * x) / ps_over_p0

            second = (1 + eta) / (1 + eta0) * Ey / chiE0 / ps_over_p0 + Bx / chiM0 - a * Bs / chiM0 / ps_over_p0
            du[4] = (1 + h * x) * second

            # du[5] = 0.0  #((1 + h*x) * (1+eta)/(1+eta0)/ps_over_p0 - 1) *
            # (-(1+eta0)/(2+eta0))

            K0 = eta0 * m_p * c^2
            # du[6] = ifelse(p.is_losing_energy, -p.energy_loss_per_meter / K0, 0.0)

            if p.is_losing_energy
                du[6] = -p.energy_loss_per_meter / K0
            else
                du[6] = 0
            end
            gamma = 1 + eta

            ps = ps_over_p0 * p_magic
            dtdS = (1 + h * x) * m_p * (1 + eta) / ps

            du[5] = dtdS
            
            px = a * p_magic 
            py = b * p_magic 

            βx = px / gamma / m_p / c
            βy = py / gamma / m_p / c
            βs = ps / gamma / m_p / c

            C₁ = (G_p + 1 / gamma)
            C₂ = G_p * gamma / (gamma + 1)
            C₃ = (G_p + (1 / (gamma + 1))) / c
            βB = (βx * Bx + βy * By + βs * Bs)

            Ωx = C₁ * Bx - C₂ * βx * βB - C₃ * (βy * Es - βs * Ey) 
            Ωy = C₁ * By - C₂ * βy * βB - C₃ * (βs * Ex - βx * Es) 
            Ωs = C₁ * Bs - C₂ * βs * βB - C₃ * (βx * Ey - βy * Ex) 

            η = 0.0

            if p.if_EDM_on
                η += p.η
            end

            C₄ = gamma / (c * (gamma + 1))
            βE = (βx * Ex + βy * Ey + βs * Es)
        
            Ωxη = η / 2 * ((βy * Bs - βs * By) + Ex / c - C₄ * βx * βE)
            Ωyη = η / 2 * ((βs * Bx - βx * Bs) + Ey / c - C₄ * βy * βE)
            Ωsη = η / 2 * ((βx * By - βy * Bx) + Es / c - C₄ * βs * βE)

            Ωx += Ωxη 
            Ωy += Ωyη 
            Ωs += Ωsη 

            Ωx *= - e / m_p * dtdS
            Ωy *= - e / m_p * dtdS
            Ωs *= - e / m_p * dtdS

            Ωy -= - h # Account for FS system rotation

            du[7] = Ωy * u[9] - Ωs * u[8]
            du[8] = Ωs * u[7] - Ωx * u[9]
            du[9] = Ωx * u[8] - Ωy * u[7]

            return nothing
        end
    end
end
