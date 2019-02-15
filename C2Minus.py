import networkx as nx
from sage.all import *
from Braid import *


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


        fully_singular = self.braid.singular_resolution()
        variables = self.R.gens()

        # Construct d1.
        # d1[I][J] will be where 1 in R/(N + LI) in cube[I] maps to in cube[J].
        d1 = {I: {} for I in self.cube.keys()}
        for i in range(len(word)):
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

    # The edge assignment for the edge (I, J) of the cube.
    def edge_assigment(self, I, J):
        answer = 0
        for i in range(len(self.word)):
            bitmask = (1 << i)
            if I & bitmask != J & bitmask:
                break
            answer += 1 if (I & bitmask) else 0
        return answer % 2




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
                    if self.graph.vdict[v][j] == None:
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
                    if x != None:
                        out_prod *= variables[x]
                for x in self.graph.vdict[v][2:]:
                    if x != None:
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
