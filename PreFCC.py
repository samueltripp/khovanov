from FCC import *
from MonomialsOfDegree import *


# represents a filtered chain complex over a polynomial ring, with chain groups built from sums of quotient rings
class PreFCC:
    def __init__(self):
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

        # would like to replace this with memoization
        genlists = {}

        # Turn vertices into sets of vertices
        i=0 # will remove after debugging
        for key, f_level in self.vertices.items():
            i+=1
            print('Processing vertex '+str(i)+'/'+str(len(self.vertices.keys())))
            if self.rings[key] not in genlists:
                print('Finding genlist')
                glist = genlist(k+1, self.rings[key])
                genlists[self.rings[key]] = [gen for grading in glist for gen in grading]
            gens = genlists[self.rings[key]]

            for g in gens:
                new_vertices[(key,g)] = f_level

        # Turn edges into sets of edges
        i=0
        for e in self.edges:
            i+=1
            print('Processing edge '+str(i)+'/'+str(len(self.edges)))
            R = self.rings[e.source]
            S = self.rings[e.target]
            S_ideal = S.ideal(0)
            if is_QuotientRing(S):
                S_ideal = S.defining_ideal()

            for x in genlists[R]:
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
