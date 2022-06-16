function get_rate_oscillatory_fit(x, y)

    linfit_res = linearfit(x, y)
    a0, a1 = linfit_res.param
    p0 = [a0, a1, 0.0 ,0.0 ,0.0]
    osc_res = curve_fit(sin_model, x, y, p0)

    return osc_res
end

function get_maximum_rate(x, y)
    delta_t, delta_y = get_moving_average(x, y, 2)
    return (delta_y[2] - delta_y[1]) / (delta_t[2] - delta_t[1])
end

function get_simple_linear_rate(x, y)
    Δx, Δy = get_moving_average(x, y, 3)
    rate1 = (Δy[2] - Δy[1]) / (Δx[2] - Δx[1])
    rate2 = (Δy[3] - Δy[2]) / (Δx[3] - Δx[2])


    return get_maximum_rate(x, y) ± abs.(rate1 - rate2)/2
end

function get_hor_precession_rate(sol)
    return get_simple_linear_rate(sol[time_index,:], sol[7,:])
end

function get_vert_precession_rate(sol)
    return get_simple_linear_rate(sol[time_index,:], sol[8,:])
end


function get_radially_compensated_vertical_precession_data(sol, sol_rad)
    # dS_y/dt = η * S_s + δ * S_x
    # dS_x/dt = a S_s
    a = get_simple_linear_rate(sol[time_index,:], sol[7,:]).val
    δ = get_simple_linear_rate(sol_rad[time_index,:], sol_rad[8,:]).val

    return sol[time_index,:], @. sol[8,:] - a * δ * sol[time_index,:]^2 / 2
end

function get_radially_compensated_4_total_vertical_precession_data(sols, sols_rad)
    t = @view sols[1][time_index,:]
    a1 = @view get_radially_compensated_vertical_precession_data(sols[1], sols_rad[1])[2][:]
    a2 = @view get_radially_compensated_vertical_precession_data(sols[2], sols_rad[2])[2][:]
    a3 = @view get_radially_compensated_vertical_precession_data(sols[3], sols_rad[3])[2][:]
    a4 = @view get_radially_compensated_vertical_precession_data(sols[4], sols_rad[4])[2][:]

    return t, @. (a1 + a2 - a3 - a4) /4
end

function get_radially_compensated_2_total_vertical_precession_data(sols, sols_rad)
    t = @view sols[1][time_index,:]
    a1 = @view get_radially_compensated_vertical_precession_data(sols[1], sols_rad[1])[2][:]
    a2 = @view get_radially_compensated_vertical_precession_data(sols[2], sols_rad[2])[2][:]
    a3 = @view get_radially_compensated_vertical_precession_data(sols[3], sols_rad[3])[2][:]
    a4 = @view get_radially_compensated_vertical_precession_data(sols[4], sols_rad[4])[2][:]

    return t, @. (a1 - a3) /2
end
