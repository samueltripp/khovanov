def genlist(maxdegree, ring):
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
                print i, ideal.reduce(prod)
                monomials.append(ideal.reduce(prod))                
    return monomials
