from FCC import *
from functools import partial
from sage.rings.quotient_ring import is_QuotientRing
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors


# represents a filtered chain complex over a polynomial ring, with chain groups built from sums of quotient rings
class PreFCC:
    def __init__(self, R):
        self.R = R
        self.vertices = {}
        self.edges = set()
        self.rings = {}

    def add_vertex(self, key, ring, f_level):
        self.vertices[key] = f_level
        self.rings[key] = ring

    def add_edge(self, source, target, coefficient):
        if coefficient != 0:
            self.edges.add(PreFCC.Edge(source, target, coefficient))

    # truncates each ring at polynomial degree k+1
    # returns an FCC over the base field
    def truncate(self, k):
        new_vertices = {}
        new_edges = set()

        gen_lists = {}

        # This is ready to be parallelized, if everything wasn't so awful.
        for I in {ring_defining_ideal(ring) for ring in self.rings.values()}:
            print('Finding genlist')
            glist = gen_list(k + 1, self.R.quotient(I))
            gen_lists[I] = [gen for grading in glist for gen in grading]

        # Turn vertices into sets of vertices
        for key, f_level in self.vertices.items():
            gens = gen_lists[ring_defining_ideal(self.rings[key])]
            for g in gens:
                new_vertices[(key,g)] = f_level

        # Turn edges into sets of edges
        i=0
        for e in self.edges:
            i+=1
            print('Processing edge '+str(i)+'/'+str(len(self.edges)))
            R = self.rings[e.source]
            S = self.rings[e.target]
            S_ideal = ring_defining_ideal(S)

            for x in gen_lists[ring_defining_ideal(R)]:
                image = e.coefficient * S.retract(R.lift(x))
                image = S_ideal.reduce(image)
                for c, y in terms(image):
                    if y.lift().degree() > k+1:
                        continue
                    new_edges.add(FCC.Edge((e.source, x), (e.target, y), c))

        return FCC(new_vertices, new_edges)

    class Edge:
        def __init__(self, source, target, coefficient):
            self.source = source
            self.target = target
            self.coefficient = coefficient

        def __repr__(self):
            return str(self.source) + '--- ' + str(self.coefficient) + ' -->' + str(self.target)


def terms(p):
    output = []
    while p != 0:
        term = (p.lc(), p.lm())
        output.append(term)
        p -= term[0]*term[1]
    return output

def ring_defining_ideal(R):
    R_ideal = R.ideal(0)
    if is_QuotientRing(R):
        R_ideal = R.defining_ideal()
    return R_ideal

def gen_list(max_degree, R):
    ideal = ring_defining_ideal(R)
    gens = R.gens()
    monomials = [[R(1)]]

    for i in range(1, max_degree + 1):
        degs = WeightedIntegerVectors(i, [1] * len(gens))
        degree_ideal = R.ideal(0)
        monomials.append([])
        for part in degs:
            prod = R(1)
            for j in range(len(gens)):
                prod *= gens[j]**part[j]
            prod = ideal.reduce(prod)
            if prod != 0 and prod not in degree_ideal:
                monomials[i].append(prod)
                degree_ideal = R.ideal(monomials[i])
    return monomials
