from sortedcontainers import SortedList
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

        self.edge_priority = {}
        for i, e in enumerate(edges):
            print ('FCC: adding edge '+str(i)+'/'+str(len(edges)))
            self.add_edge(e)

        # set up the proper priorities before any reduction
        self.edges = {}
        for i, e in enumerate(edges):
            print ('FCC: updating edge priority '+str(i)+'/'+str(len(edges)))
            self.edge_priority[e] = len(self.inv[e.target]) * len(self.outv[e.source])
            if self.delta_f(e) in self.edges:
                self.edges[self.delta_f(e)].add((e, self.edge_priority[e]))
            else:
                self.edges[self.delta_f(e)] = SortedList([(e, self.edge_priority[e])], key = lambda x: x[1])

    # reduces this chain complex until there are no edges remaining
    # if a page is specified, then stop reducing at that page of the spectral sequence
    def reduce(self, page=None):
        for i in itertools.count():
            if (page is not None and i == page) or self.num_edges == 0:
                return

            while i in self.edges and len(self.edges[i]) > 0:
                u = self.edges[i].pop(index = 0)[0]
                self.num_edges -= 1

                x, y, c = u.source, u.target, u.coefficient

                del self.inv[y][x]
                del self.outv[x][y]

                new_edges = []
                print ('Reduction at '+str(i)+': '+str(len(self.vertices))+' vertices, '+
                        str(len(self.edges[i]))+' edges')
                for w in self.inv[y]:
                    t = -self.inv[y][w].coefficient * 1 / c
                    for z in self.outv[x]:
                        e = self.get_edge(w, z)
                        e.coefficient += t * self.outv[x][z].coefficient
                        new_edges.append(e)

                for e in new_edges:
                    self.add_edge(e)

                self.remove_vertex(x)
                self.remove_vertex(y)

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

        t = self.edge_priority[e] if e in self.edge_priority else None
        self.edge_priority[e] = len(self.inv[e.target]) * len(self.outv[e.source])

        if self.delta_f(e) in self.edges:
            if (e, t) in self.edges[self.delta_f(e)]:
                self.edges[self.delta_f(e)].remove((e, t))
                self.num_edges -= 1
            self.edges[self.delta_f(e)].add((e, self.edge_priority[e]))
        else:
            self.edges[self.delta_f(e)] = SortedList([(e, self.edge_priority[e])], key = lambda x: x[1])
        self.num_edges += 1

    def remove_edge(self, e):
        if e in self.edge_priority:
            if (e, self.edge_priority[e]) in self.edges[self.delta_f(e)]:
                self.edges[self.delta_f(e)].remove((e, self.edge_priority[e]))
                self.num_edges -= 1
            del self.edge_priority[e]
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
