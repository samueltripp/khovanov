from sage.rings.quotient_ring import is_QuotientRing
from sage.all import *

def genlist(maxdegree, ring):
    ideal = ring.ideal(0)
    if is_QuotientRing(ring):
        ideal = ring.defining_ideal()
    gens = ring.gens()
    monomials = [[ring(1)]]

    for i in range(1, maxdegree + 1):
        degs = WeightedIntegerVectors(i, [1] * len(gens))
        degreeideal = ring.ideal(0)
        monomials.append([])
        for part in degs:
            prod = ring(1)
            for j in range(len(gens)):
                prod *= gens[j]**part[j]
            prod = ideal.reduce(prod)
            if prod != 0 and prod not in degreeideal:
                monomials[i].append(prod)
                degreeideal = ring.ideal(monomials[i])
    return monomials
