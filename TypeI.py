class TypeI:
    # returns the graph S'_{2n} in (V,E) format
    # for now, vertices are strings with prefixes out-, top-, mid-, bot-, and in-
    # edges are unlabeled
    @staticmethod
    def graph(n):
        V = {'out'+str(i) for i in range(1, 2*n+1)} | \
            {'top'+str(i) for i in range(1, n+1)} | \
            {'mid'+str(i) for i in range(1, n)} | \
            {'bot'+str(i) for i in range(1, n)} | \
            {'in'+str(i) for i in range(1, 2*n+1)}

        E = {('top'+str((i+1)//2), 'out'+str(i)) for i in range(1, 2*n+1)} | \
            {('mid' + str((i-1)//2), 'top' + str(i//2)) for i in range(3, 2*n+1)} | \
            {('bot' + str((i-1)//2), 'mid' + str(i//2)) for i in range(3, 2*n)} | \
            {('in' + str(i+1), 'bot' + str(i//2)) for i in range(2, 2*n)} | \
            {('in1', 'top1'), ('in2', 'mid1'), ('bot'+str(n-1), 'top'+str(n))}

        return V, E
