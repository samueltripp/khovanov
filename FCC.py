

class FCC:
    def __init__(self, vertices, edges):
        self.vertices = vertices
        self.inv = {v: {} for v in vertices}
        self.outv = {v: {} for v in vertices}
        self.unit_edges = set()
        
        for e in edges:
            self.add_edge(e)

    def reduce(self):
        while len(self.unit_edges) > 0:
            u = self.unit_edges.pop()
            x, y, c = u.source, u.target, u.coefficient

            new_edges = []
            for w in self.inv[u.target]:
                for z in self.outv[u.source]:
                    if w != u.source and z != u.target:
                        e = self.get_edge(w,z)
                        e.coefficient += -self.get_edge(w, y).coefficient * c * self.get_edge(x, z).coefficient
                        new_edges.append(e)

            for e in new_edges:
                self.add_edge(e)

            self.remove_vertex(x)
            self.remove_vertex(y)

    def get_edge(self, source, target):
        if target in self.outv[source]:
            return self.outv[source][target]
        else:
            return FCC.Edge(source, target, 0)

    def add_edge(self, e):
        if e.coefficient == 0:
            self.remove_edge(e)
            return
        elif e.coefficient == 1 or e.coefficient == -1:
            self.unit_edges.add(e)
        else:
            if e in self.unit_edges:
                self.unit_edges.remove(e)

        self.inv[e.target][e.source] = e
        self.outv[e.source][e.target] = e

    def remove_edge(self, e):
        if e in self.unit_edges:
            self.unit_edges.remove(e)
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

        self.vertices.remove(v)

    class Edge:
        def __init__(self, source, target, coefficient):
            self.source = source
            self.target = target
            self.coefficient = coefficient

        def __repr__(self):
            return str(self.source) + '--' + str(self.coefficient) + '--' + str(self.target)
