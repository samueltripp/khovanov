from sage.all import *
import networkx as nx


class C2Minus:
    # resolution is a Braid.Graph object.
    def __init__(self, resolution):
        self.graph = resolution
        # Base ring.
        self.R = PolynomialRing(QQ, ['U%s'%p for p in range(1, 
            len(self.graph.edges) + 1)])

        print self.N()
        print self.LI()

        self.complex = self.R.quotient(self.N() + self.LI())
        # Tensor this with the chain complexes, L_D^+. 

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
