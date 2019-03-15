from FCC import *
from functools import partial
from multiprocessing import Pool
from sage.rings.quotient_ring import is_QuotientRing
from sage.combinat.integer_vector_weighted import WeightedIntegerVectors
import random


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
        # for i, I in enumerate({ring_defining_ideal(ring) for ring in self.rings.values()}):
        #     glist = gen_list(k + 1, self.R.quotient(I), str(i))
        #     gen_lists[I] = [gen for grading in glist for gen in grading]

        pool = Pool()
        ideal_list = list({ring_defining_ideal(ring) for ring in self.rings.values()})
        quotient_rings = [self.R.quotient(ideal) for ideal in ideal_list]
        for (idx, lifted_gens) in pool.imap_unordered(partial(gen_list_helper, k), enumerate(quotient_rings)):
            gen_lists[ideal_list[idx]] = [quotient_rings[idx].retract(gen) for gen in lifted_gens]

        pool.close()
        pool.join()

        # Turn vertices into sets of vertices
        for key, f_level in self.vertices.items():
            gens = gen_lists[ring_defining_ideal(self.rings[key])]
            for g in gens:
                new_vertices[(key,g)] = f_level

        # Turn edges into sets of edges
        i=0
        for e in self.edges:
            i+=1
            print('PreFCC: processing edge '+str(i)+'/'+str(len(self.edges)))
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

def gen_list(max_degree, R, ideal_idx):
    ideal = ring_defining_ideal(R)
    gens = R.ambient().gens()
    monomials = [[R(1)]]

    for i in range(1, max_degree + 1):
        degs = WeightedIntegerVectors(i, [1] * len(gens))
        degree_ideal = R.ideal(0)
        monomials.append([])
        for m, part in enumerate(degs):
            print('gen_list for '+str(ideal_idx)+': degree '+str(i)+'/'+str(max_degree)+
                ', index '+str(m)+'/'+str(len(degs)))
            prod = R.ambient()(1)
            for j in range(len(gens)):
                prod *= gens[j]**part[j]
            p2 = ideal.reduce(R.retract(prod))

            if prod not in degree_ideal:
                monomials[i].append(prod)
                degree_ideal += R.ideal(prod)

    return monomials

def gen_list_piece(i,R,ideal_idx,ideal,gens):
    output = []
    partition_list = WeightedIntegerVectors(i,[1]*len(gens))
    degs = []
    for part in partition_list:
        degs.append(part)

    degs.reverse()
    degree_ideal = R.ideal(0)

    for m, part in enumerate(degs):
        print('gen_list for '+str(ideal_idx)+': degree '+str(i)+', index '+str(\
m)+'/'+str(len(degs)))
        prod = R.ambient()(1)
        for j in range(len(gens)):

        p2 = ideal.reduce(R.retract(prod))
        print p2
        test = degree_ideal.reduce(p2).is_zero()
        if not test:
            output.append(p2)
            degree_ideal += R.ideal(p2)
            print p2

    return output



def gen_list_helper(k, (ideal_idx, R)):
    output = [[]]
    ideal = ring_defining_ideal(R)
    gens = R.ambient().gens()
    for i in range(1,k+2):
        output.append(gen_list_piece(i,R,ideal_idx,ideal,gens))
    glist = output
    #glist = gen_list(k + 1, R, ideal_idx)                                                        
    # Return the .lift() because for some reason, sage bugs out if you return                     
    # anything that has to do with the defining ideal of R.                                       
    return (ideal_idx, [gen.lift() for grading in glist for gen in grading])
