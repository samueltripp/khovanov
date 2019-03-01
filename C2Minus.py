import networkx as nx
from sage.all import *
from Braid import *
from PreFCC import *


# C_2^- for a braid.
class C2Minus:
    # Same as in Braid:
    # n - *half* the number of strands in the braid
    # word - a word in the braid group B_{2n}, given as a tuple {-(2n-1),...,-1,1,...,(2n-1)}*
    def __init__(self, n, word):
        self.n = n
        self.word = word
        self.braid = Braid(n, word)

        self.R = PolynomialRing(QQ, ['U%s' % (p + 1) for p in range(6*n-4 + 2*len(word))])

        self.cube = {key: C2Minus.ResolutionComplex(self.R, v) for
                            key, v in self.braid.cube_of_resolutions().items()}

        self.LDPlus = self.initLDPlus()

        self.d1 = self.construct_d1()

    # Construct d1.
    # d1[I][J] will be where 1 in R/(N + LI) in cube[I] maps to in cube[J].
    def construct_d1(self):
        fully_singular = self.braid.singular_resolution()
        variables = self.R.gens()

        d1 = {I: {} for I in self.cube.keys()}
        for i in range(len(self.word)):
            bitmask = (1 << i)

            mult = 1
            if self.word[i] < 0:
                # For a negative crossing, need to multiply by a factor.
                edges = fully_singular.vdict['b' + str(i + 1)]
                mult = variables[edges[1]] - variables[edges[3]]

            for I, v in self.cube.items():
                # Look at all the vertices where we can increase height by flipping i.
                if (I & bitmask): continue
                J = I | bitmask

                # TODO: tensor mult with the identity when we have an actual chain complex.
                d1[I][J] = ((-1)**self.edge_assigment(I, J)) * self.cube[J].complex(mult)

        return d1

    # The edge assignment for the edge (I, J) of the cube.
    def edge_assigment(self, I, J):
        answer = 0
        for i in range(len(self.word)):
            bitmask = (1 << i)
            if I & bitmask != J & bitmask:
                break
            answer += 1 if (I & bitmask) else 0
        return answer % 2

    def initLDPlus(self):
        variables = self.R.gens()
        maps = []
        vdict = self.braid.singular_resolution().vdict

        for v in vdict.keys():
            if (v[:2] == 'to') or (v[:2] == 'mi') or (v[:2] == 'bo'):
                insum, outsum = 0, 0
                for x in vdict[v][:2]:
                    if x is not None:
                        outsum += variables[x]
                for x in vdict[v][2:]:
                    if x is not None:
                        insum += variables[x]
                maps.append(outsum - insum)
                maps.append(outsum + insum)
                
        matrix_maps = constructmatrix(maps,self.R)
        return matrix_maps

    # returns C2Minus as a PreFCC
    # vertex keys are pairs (crossing_key, s_key)
    def preFCC(self):
        pfcc = PreFCC()
        # the size of the cubes at each vertex of the crossing cube
        m = self.LDPlus[0].nrows()

        for crossing_key in self.cube:
            Q = self.R.quotient(
                    self.cube[crossing_key].N() + self.cube[crossing_key].LI() + (self.R.gen(0) * self.R.gen(0),))
            for s_key in range(m):
                for z2_grading in range(2):
                    crossing_height = bin(crossing_key).count("1")

                    # add vertices
                    pfcc.add_vertex((crossing_key, s_key, z2_grading), Q, crossing_height)

                    # add d0 edges
                    for target in range(m):
                        pfcc.add_edge(
                            (crossing_key, s_key, z2_grading),
                            (crossing_key, target, 1 - z2_grading),
                            self.LDPlus[z2_grading][s_key, target])

                    # add d1 edges
                    for target, coefficient in self.d1[crossing_key].items():
                        pfcc.add_edge(
                            (crossing_key, s_key, z2_grading),
                            (target, s_key, z2_grading),
                            coefficient)

        return pfcc

    # returns the homology in polynomial degree <= k
    def homology(self, k):
        print('Computing preFCC...')
        pfcc = self.preFCC()
        print('Computing truncated FCC...')
        fcc = pfcc.truncate(k)
        print('Reducing to E^1 page...')
        fcc.reduce(page=1)
        print('Truncating E^1 page...')
        fcc.truncate(k)
        print('Reducing to E^infty page...')
        fcc.reduce()
        print('Done!')
        return fcc

    # C_2^-(S) for a complete resolution S.
    class ResolutionComplex:
        # R - the base polynomial ring.
        # resolution - a Braid.Graph object.
        def __init__(self, R, resolution):
            self.R = R
            self.graph = resolution
            self.complex = self.R.quotient(self.N() + self.LI())
            # TODO: tensor this with the chain complexes, L_D^+.

        # Given a cycle as a list of nodes, find all possible edge sets that form
        # that cycle. This is needed because might have multiple edges between two
        # nodes, and nx doesn't handle that.
        def edge_resolutions(self, cycle):
            if len(cycle) <= 1:
                return [[]]

            result = []

            for e in self.graph.vdict[cycle[0]][:2]:
                if e != None and self.graph.edges[e][1] == cycle[1]:
                    # This edge works for the first two nodes.
                    result.extend([[e] + x for x in self.edge_resolutions(cycle[1:])])

            return result

        # Compute the \prod_{u going out} u - \prod_{u going in} u for the cycle.
        def cycle_relation(self, cycle):
            variables = self.R.gens()

            in_prod, out_prod = 1, 1
            for i in range(len(cycle) - 1):
                e0, e1 = cycle[i], cycle[i + 1]
                # Where the edge points to.
                v = self.graph.edges[e0][1]

                # Find the index of the incoming edge.
                j = 0
                while self.graph.vdict[v][j] != e0:
                    j += 1

                # Go over every edge "outside" the cycle.
                k = len(self.graph.vdict[v])
                while self.graph.vdict[v][(j + 1) % k] != e1:
                    j = (j + 1) % k
                    if self.graph.vdict[v][j] is None:
                        continue

                    # The first half is things going out.
                    if j < 2:
                        out_prod *= variables[self.graph.vdict[v][j]]
                    else:
                        in_prod *= variables[self.graph.vdict[v][j]]

            return out_prod - in_prod

        # Compute the ideal N for the resolution.
        def N(self):
            variables = self.R.gens()
            ideal_generators = []

            # Get the quadratic relations from the vertices.
            for v in self.graph.vdict.keys():
                in_prod, out_prod = 1, 1
                for x in self.graph.vdict[v][:2]:
                    if x is not None:
                        out_prod *= variables[x]
                for x in self.graph.vdict[v][2:]:
                    if x is not None:
                        in_prod *= variables[x]
                ideal_generators.append(out_prod - in_prod)

            # Use networkx for getting all cycles.
            G = nx.DiGraph()
            # Don't add the first edge.
            G.add_edges_from(self.graph.edges[1:])
            for cycle in nx.simple_cycles(G):
                cycle.append(cycle[0])
                for edge_cycle in self.edge_resolutions(cycle):
                    edge_cycle.append(edge_cycle[0])
                    ideal_generators.append(self.cycle_relation(edge_cycle))

            return self.R.ideal(ideal_generators)

        # The LI bit.
        def LI(self):
            variables = self.R.gens()
            ideal_generators = []

            for v in self.graph.vdict.keys():
                if v[0] == 'b' and v[1].isdigit():
                    relation = 0
                    for x in self.graph.vdict[v][:2]:
                        if x != None:
                            relation += variables[x]
                    for x in self.graph.vdict[v][2:]:
                        if x != None:
                            relation -= variables[x]
                    ideal_generators.append(relation)

            return self.R.ideal(ideal_generators)

            
def constructmatrix(maps,R):
    if len(maps) == 2:
        L = matrix(R,[[maps[0]]])
        Lplus = matrix(R,[[maps[1]]])
        return [L,Lplus]
    if len(maps) == 4:
        L = matrix(R,[[maps[2],maps[1]],[maps[0],maps[3]]])
        Lplus = matrix(R, [[maps[3],maps[1]],[maps[0],maps[2]]])
        return [L,Lplus]
    else:
        nummaps = len(maps)
        D1, D2 = constructmatrix(maps[:nummaps-2],R)
        diag1 = maps[nummaps-2]*identity_matrix(R,D1.nrows())
        diag2 = maps[nummaps-1]*identity_matrix(R,D1.nrows())
        L = block_matrix(R,[[diag1,D2],[D1,diag2]])
        LPlus = block_matrix(R, [[diag2,D2],[D1,diag1]])
        return [L,LPlus]
