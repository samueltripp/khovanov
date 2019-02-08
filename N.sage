import networkx as nx

def edge_resolutions(cycle, graph):
    # Given a cycle as a list of nodes, find all possible edge sets that form
    # that cycle. This is needed because might have multiple edges between two
    # nodes.

    if len(cycle) <= 1:
        return [[]]

    result = []

    edges, inv, outv, _ = graph
    for e in outv[cycle[0]]:
        if edges[e][1] == cycle[1]:
            # This edge works.
            result.extend([[e] + x for x in edge_resolutions(cycle[1:], graph)])

    return result

def cycle_relation(R, graph, cycle):
    variables = R.gens()
    edges, inv, outv, allv = graph

    in_prod, out_prod = 1, 1
    for i in range(len(cycle) - 1):
        e0, e1 = cycle[i], cycle[i + 1]
        # Where the edge points to.
        v = edges[e0][1]

        # Find the index of the incoming edge.
        j = 0
        while allv[v][j] != e0:
            j += 1

        # Go over every edge "outside" the cycle.
        k = len(allv[v])
        while allv[v][(j + 1) % k] != e1:
            j = (j + 1) % k
            # Bad shikhin?
            if j < 2:
                out_prod *= variables[allv[v][j]]
            else:
                in_prod *= variables[allv[v][j]]

    return out_prod - in_prod

def N(R, graph):
    variables = R.gens()
    edges, inv, outv, _ = graph
    ideal_generators = []

    # Get the quadratic relations from the vertices.
    for v in inv.keys():
        in_prod, out_prod = 1, 1
        for x in inv[v]:
            in_prod *= variables[x]
        for x in outv[v]:
            out_prod *= variables[x]
        ideal_generators.append(out_prod - in_prod)

    # Use networkx for getting all cycles.
    G = nx.DiGraph()
    G.add_edges_from(edges)
    for cycle in nx.simple_cycles(G):
        cycle.append(cycle[0])
        for edge_cycle in edge_resolutions(cycle, graph):
            edge_cycle.append(edge_cycle[0])
            ideal_generators.append(cycle_relation(R, graph, edge_cycle))

    return ideal_generators

vertices = list(range(4))
edges = [(0, 0), (0, 1), (1, 0), (2, 3), (2, 3), (1, 2), (3, 1), (3, 2)]

inv = {v: [] for v in vertices}
outv = {v: [] for v in vertices}
for i, edge in enumerate(edges):
    outv[edge[0]].append(i)
    inv[edge[1]].append(i)
inv = {v: o[::-1] for v, o in inv.items()}
allv = {v: outv[v] + inv[v] for v in inv.keys()}

graph = (edges, inv, outv, allv)
R = PolynomialRing(QQ, 'u', len(edges))

N_my = R.ideal(N(R, graph))

R.inject_variables()
N_nate = R.ideal([u2 - u1, u5 - u6, u3 * u4 - u5 * u7])

print "The real ideal:", N_nate
print "My presentation:", N_my
print "Equal?", N_nate == N_my
