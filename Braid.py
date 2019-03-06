class Braid:
    # n - *half* the number of strands in the braid
    # word - a word in the braid group B_{2n}, given as a tuple {-(2n-1),...,-1,1,...,(2n-1)}*
    #
    # conventions: braids are pictured vertically, strands are numbered 1-2n left-to-right,
    # and the generator i (-i) moves the ith strand under (over) the (i+1)st strand.
    def __init__(self, n, word):
        for i in word:
            assert 1 <= abs(i) < 2*n, "Invalid word"

        self.n = n
        self.word = word

    # returns a dict representing image of this braid under the canonical projection B_{2n} -> S_{2n}
    def permutation(self):
        p = {i: i for i in range(1, 2*self.n+1)}
        for swap in reversed(self.word):
            swap = abs(swap)
            p[swap], p[swap+1] = p[swap+1], p[swap]
        return p

    # returns the fully singular resolution of S'_{2n}+B
    # output is a Braid.Graph object
    # for now, vertices are strings with prefixes top-, mid-, bot-, b-, br-, and bl-
    def singular_resolution(self):
        n = self.n
        word = self.word

        g = Braid.Graph(6*n-4 + 2*len(word))

        if n >= 2:
            for i in range(1, n):
                g.add_edge(2*n+2*i-1, 'mid'+str(i), 1, 'top'+str(i), 3)
                g.add_edge(2*n+2*i, 'mid'+str(i), 2, 'top'+str(i+1), 4)

            for i in range(1, n-1):
                g.add_edge(4*n+2*i-3, 'bot'+str(i), 1, 'mid'+str(i), 3)
                g.add_edge(4*n+2*i-2, 'bot'+str(i), 2, 'mid'+str(i+1), 4)

            g.add_edge(6*n-5, 'bot'+str(n-1), 1, 'mid'+str(n-1), 3)
            g.add_edge(6*n-4, 'bot'+str(n-1), 2, 'top'+str(n), 3)

            strands = [None, ('top1', 4), ('mid1', 4)]
            for i in range(1, n):
                strands.append(('bot'+str(i), 4))
                strands.append(('bot'+str(i), 3))
            next_strand_number = 6*n-3
        else:
            strands = [None, ('top1', 4), ('top1', 3)]
            next_strand_number = 3

        for i, swap in enumerate(word):
            i=i+1
            s = abs(swap)

            g.add_edge(next_strand_number, 'b'+str(i), 1, *strands[s])
            g.add_edge(next_strand_number+1, 'b'+str(i), 2, *strands[s+1])
            next_strand_number += 2
            strands[s] = ('b'+str(i), 4)
            strands[s+1] = ('b'+str(i), 3)

        for i in range(1, n+1):
            g.add_edge(2*i-1, 'top'+str(i), 1, *strands[2*i-1])
            g.add_edge(2*i, 'top'+str(i), 2, *strands[2*i])

        return g

    # returns a dictionary with integer keys (thought of as binary strings) and resolutions (type Braid.Graph) as values
    # time: O(c * 2^c)
    def cube_of_resolutions(self):
        k = len(self.word)
        cube = {}

        initial_key = int(''.join(['0' if s > 0 else '1' for s in reversed(self.word)]), 2)
        initial_resolution = self.singular_resolution()

        frontier = {initial_key: initial_resolution}
        new = {}

        for _ in range(0,k):
            for key, resolution in frontier.items():
                for i in range(0,k):
                    if key ^ 2**i in new or 'b'+str(i+1) not in resolution.vdict:
                        continue
                    else:
                        new[key ^ 2**i] = resolution.split_vertex(i+1)
            cube.update(frontier)
            frontier = new
        cube.update(frontier)

        return cube

    # represents the enhanced graph structure needed to represent the braid diagrams
    # attributes:
    #   edges - a list of edges
    #   vdict - a dictionary {v:(e1,e2,e3,e4)} where the e_i are edge indices, ordered clockwise around the vertex with
    #           e1 and e2 leaving, e3 and e4 entering
    class Graph:
        def __init__(self, num_edges):
            self.edges = [None]*num_edges
            self.vdict = {}

        def add_edge(self, var, source, source_position, target, target_position):
            e = (source, target)

            assert self.edges[var-1] is None, "Edge variable is already taken"
            self.edges[var-1] = e

            if source not in self.vdict:
                self.vdict[source] = [None, None, None, None]
            else:
                assert self.vdict[source][source_position-1] is None, 'Source position is already occupied'
            self.vdict[source][source_position-1] = var-1

            if target not in self.vdict:
                self.vdict[target] = [None, None, None, None]
            else:
                assert self.vdict[target][target_position - 1] is None, 'Target position is already occupied'
            self.vdict[target][target_position - 1] = var-1

        # returns the graph given by replacing braid vertex b_i with the oriented smoothing,
        # leaving behind two vertices bl_i and br_i
        # time: O(c)
        def split_vertex(self, i):
            edges = list(self.edges)
            vdict = dict(self.vdict)

            v = 'b'+str(i)
            assert v in vdict, "Vertex does not exist"

            vl = 'bl'+v[1:]
            vr = 'br'+v[1:]

            i1 = vdict[v][0]
            e1 = edges[i1]
            i2 = vdict[v][1]
            e2 = edges[i2]
            i3 = vdict[v][2]
            e3 = edges[i3]
            i4 = vdict[v][3]
            e4 = edges[i4]
            edges[i1] = (vl, e1[1])
            edges[i2] = (vr, e2[1])
            edges[i3] = (e3[0], vr)
            edges[i4] = (e4[0], vl)

            vdict[vl] = [i1, None, None, i4]
            vdict[vr] = [None, i2, i3, None]
            del vdict[v]

            g = Braid.Graph(0)
            g.edges = edges
            g.vdict = vdict
            return g
