function get_solution_args(u₀, p; cw=true, ode=bmt_equation, kwargs...)
    u = deepcopy(u₀)

    if !cw
        u[3] = u[3] * (-1)
        u[4] = u[4] * (-1)
        u[8] = u[8] * (-1)
    end
    p = deepcopy(p)
    ending_time = p.turns * p.ring.total_length
    time_span = (p.starting_time, ending_time)

    ode_func = if cw
        (du, u, p, t) -> ode(du, u, p, t, cw=true)
    else
        (du, u, p, t) -> ode(du, u, p, t, cw=false)
    end

    prob = ODEProblem(ode_func, u, time_span, p)

    changeRegionsCallback = IterativeCallback((integrator) -> next_region_time(integrator),
                                              (integrator) -> p.region_change_function!(integrator, cw),
                                              save_positions = (p.save_positions,p.save_positions),)

    RFKickingCallback = IterativeCallback(
        (integrator) -> next_rf_time(integrator, cw),
        (integrator) -> apply_rf!(integrator, cw),
        save_positions = (p.save_positions,p.save_positions),)

    newkwargs = Dict(
        :callback => p.RF_on ? CallbackSet(changeRegionsCallback, RFKickingCallback) : changeRegionsCallback,
        :dense    => false,
        :saveat   => p.saving_time_step,
        :reltol   => p.state_tolerance_rel,
        :abstol   => p.state_tolerance_abs,
        :dtmax    => p.internal_max_time_step,
        :maxiters => 1e10,
        kwargs...)

    args = prob, p.alg

    return args, newkwargs
end

function get_solution(u₀, p; cw=true, kwargs...)

    args, newkwargs = get_solution_args(u₀, p; cw=cw, kwargs...)
    sol = solve(args...; newkwargs...)


    if !cw
        sol[8,:] .= -sol[8,:]
        sol[3,:] .= -sol[3,:]
        sol[4,:] .= -sol[4,:]
    end


    return sol
end
