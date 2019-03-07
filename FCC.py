import itertools


# represents a filtered chain complex over a field
class FCC:
    # vertices - a dictionary {v : filtration level}
    # edges - a set [FCC.Edge]
    def __init__(self, vertices, edges):
        self.vertices = vertices
        self.inv = {v: {} for v in vertices}
        self.outv = {v: {} for v in vertices}
        self.edges = {}  # a dictionary {filtration level change : {FCC.Edge}}
        self.num_edges = 0  # keeps track of the number of edges in self.edges

        for e in edges:
            self.add_edge(e)

    # reduces this chain complex until there are no edges remaining
    # if a page is specified, then stop reducing at that page of the spectral sequence
    def reduce(self, page=None):
        for i in itertools.count():
            if (page is not None and i == page) or self.num_edges == 0:
                return

            while i in self.edges and len(self.edges[i]) > 0:
                # j = self.find_good_edge(self,i)
                u = self.edges[i].pop()
                self.num_edges -= 1
                x, y, c = u.source, u.target, u.coefficient

                new_edges = []
                for w in self.inv[u.target]:
                    if w == u.source:
                        continue

                    t = -self.get_edge(w, y).coefficient * 1 / c
                    for z in self.outv[u.source]:
                        if z != u.target:
                            e = self.get_edge(w, z)
                            e.coefficient += t * self.get_edge(x, z).coefficient
                            new_edges.append(e)

                for e in new_edges:
                    self.add_edge(e)

                self.remove_vertex(x)
                self.remove_vertex(y)

    # find an edge that's connected to a vertex with small degree
    def find_good_edge(self,i):
        v_degs = {e.source: len(inv(e.source))+len(outv(e.source)) for e in self.edges[i]}
        return min(v_degs,key = v_degs.get)

    # remove all vertices with polynomial degree > k
    def truncate(self, k):
        for v in self.vertices.keys():
            if v[1].lift().degree() > k:
                self.remove_vertex(v)

    def get_edge(self, source, target):
        if target in self.outv[source]:
            return self.outv[source][target]
        else:
            return FCC.Edge(source, target, 0)

    def add_edge(self, e):
        if e.coefficient == 0:
            self.remove_edge(e)
            return

        self.inv[e.target][e.source] = e
        self.outv[e.source][e.target] = e
        if self.delta_f(e) in self.edges:
            if e in self.edges[self.delta_f(e)]:
                self.edges[self.delta_f(e)].remove(e)
                self.num_edges -= 1
            self.edges[self.delta_f(e)].add(e)
        else:
            self.edges[self.delta_f(e)] = {e}
        self.num_edges += 1

    def remove_edge(self, e):
        if e in self.edges[self.delta_f(e)]:
            self.edges[self.delta_f(e)].remove(e)
            self.num_edges -= 1
        if e.source in self.inv[e.target]:
            del self.inv[e.target][e.source]
        if e.target in self.outv[e.source]:
            del self.outv[e.source][e.target]

    def remove_vertex(self, v):
        for e in self.inv[v].values():
            self.remove_edge(e)
        for e in self.outv[v].values():
            self.remove_edge(e)

        del self.inv[v]
        del self.outv[v]

        del self.vertices[v]

    # returns how much this edge increases the filtration level
    def delta_f(self, e):
        return self.vertices[e.target] - self.vertices[e.source]

    class Edge:
        def __init__(self, source, target, coefficient):
            self.source = source
            self.target = target
            self.coefficient = coefficient

        def __repr__(self):
            return str(self.source) + '--- ' + str(self.coefficient) + ' -->' + str(self.target)

        def __eq__(self, other):
            return self.source == other.source and self.target == other.target

        def __hash__(self):
            return hash((self.source, self.target))
