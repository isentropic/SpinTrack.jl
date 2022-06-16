function reverse_polarity(p)
    p = deepcopy(p)
    quads = filter(i -> i isa MagneticQuadrupole, p.ring.ringElements)
    for q in quads
        q.k = -q.k
    end
    return p
end
