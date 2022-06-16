function get_total_fields(u, p, t, ringElement; cw=true)::Tuple{StaticArrays.SVector{3, Float64}, StaticArrays.SVector{3, Float64}, Float64}
    E, B, V = getFields(u, p, t, ringElement)
    # return zeros(SVector{3}), zeros(SVector{3}), 0.0
    # return E,B,V
    Eg, Bg, Vg = get_global_fields(u, p, t)

    return E+Eg, B+Bg, V+Vg
end

function get_total_fields!(ebv, u, p, t, ringElement; cw=true)
    fill!(ebv, 0)
    get_global_fields!(ebv, u, p, t)
    get_fields!(ebv, u, p, t, ringElement)
    return nothing
end

function get_global_fields(u, p, t)
    x, y = u[1], u[3]

    Efield = zeros(SVector{3})
    Bfield = zeros(SVector{3})
    V = 0.0

    if p.global_E_y
        Efield += SA[0, p.E_y, 0]
        V += p.E_y * y
    end

    if p.global_B_R
        # B_R Multipoles
        B_R = p.B_R
        if p.B_multipole != 0
            s = t % p.ring.total_length
            fractional_position = s/p.ring.total_length
            B_R = B_R * cos(fractional_position * p.B_multipole * 2 * pi + p.B_R_phase)
        end
        Bfield += SA[B_R, 0, 0]
        # Bfield = setindex(Bfield, Bfield[1] + B_R, 1)
    end

    if p.global_B_y
        Bfield += SA[0, p.B_y, 0]
    end

    if p.global_B_L
        Bfield += SA[0, 0, p.B_L]
    end


    return Efield, Bfield, V
end

function get_global_fields!(ebv, u, p, t)
    x, y = u[1], u[3]

    if p.global_E_y
        ebv[2] += p.E_y
        ebv[7] += p.E_y * y
    end

    if p.global_B_R
        # B_R Multipoles
        B_R = p.B_R
        if p.B_multipole != 0
            s = t % p.ring.total_length
            fractional_position = s/p.ring.total_length
            B_R = B_R * cos(fractional_position * p.B_multipole * 2 * pi + p.B_R_phase)
        end
        ebv[4] += B_R
        # Bfield = setindex(Bfield, Bfield[1] + B_R, 1)
    end

    if p.global_B_y
        ebv[5] += p.B_y
    end

    if p.global_B_L
        ebv[6] += p.B_L
    end
    return nothing
end
