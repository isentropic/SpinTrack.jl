function reverse_polarity(p)
    p = deepcopy(p)
    quads = filter(i -> i isa MagneticQuadrupole, p.ring.ringElements)
    for q in quads
        q.k = -q.k
    end
    return p
end

function get_4_solutions(u, p, kwargs...)
    p2 = reverse_polarity(p)
    # println("Turns: $(p.noOfpredictedTurns), all 4 cases")
    sols = Array{ODESolution}(undef, 4)
    args = ((u,p), (u,p2), (u, p), (u, p2))
    iscw = (true, true, false, false)
    for i = 1:4
        sols[i] = get_solution(args[i]..., cw=iscw[i], kwargs...)
    end

    return sols
end

function get_4_solutions_threaded(u, p, kwargs...)
    p2 = reverse_polarity(p)
    # println("Turns: $(p.noOfpredictedTurns), all 4 cases")
    sols = Array{ODESolution}(undef, 4)
    args = ((u,p), (u,p2), (u, p), (u, p2))
    iscw = (true, true, false, false)
    Threads.@threads for i = 1:4
        sols[i] = get_solution(args[i]..., cw=iscw[i], kwargs...)
    end

    return sols
end

