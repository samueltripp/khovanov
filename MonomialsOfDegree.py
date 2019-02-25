from sage.rings.quotient_ring import is_QuotientRing

def genlist(maxdegree, ring):
    ideal = ring.ideal(0)
    if is_QuotientRing(ring):
        ideal = ring.defining_ideal()
    gens = ring.gens()
    monomials = [ring(1)]

    for i in range(maxdegree+1):
        degs = WeightedIntegerVectors(i,[1]*len(gens))
        for part in degs:
            prod = ring(1)
            for j in range(len(gens)):
                prod *= gens[j]**part[j]
            if ideal.reduce(prod) != 0 and ideal.reduce(prod) not in monomials:
                monomials.append(ideal.reduce(prod))                
    return monomials
