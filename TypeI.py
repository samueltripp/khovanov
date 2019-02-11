class TypeI:
    # returns the graph S'_{2n}+B in (V,E) format
    # for now, vertices are strings with prefixes out-, top-, mid-, bot-, and in-
    # edges are unlabeled, but ordered such that E[i] corresponds to U_i
    @staticmethod
    def graph(braid):
        n = braid.n//2
        p = braid.permutation()

        V = {'top'+str(i) for i in range(1, n+1)} | \
            {'mid'+str(i) for i in range(1, n)} | \
            {'bot'+str(i) for i in range(1, n)}

        E = []
        for i in range(1, 2*n+1):
            if p[i] == 1:
                E.append(('top'+str((i+1)//2), 'top1'))
            elif p[i] == 2:
                E.append(('top'+str((i+1)//2), 'mid1'))
            else:
                E.append(('top'+str((i+1)//2), 'bot'+str((p[i]-1)//2)))
        for i in range(1, n):
            E.append(('mid'+str(i), 'top'+str(i)))
            E.append(('mid'+str(i), 'top'+str(i+1)))
        for i in range(1, n-1):
            E.append(('bot'+str(i), 'mid'+str(i)))
            E.append(('bot'+str(i), 'mid'+str(i+1)))
        E.append(('bot'+str(n-1), 'mid'+str(n-1)))
        E.append(('bot'+str(n-1), 'top'+str(n)))

        return V, E
