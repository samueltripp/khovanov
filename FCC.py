import numpy as np


class FCC:
    def __init__(self, m):
        self.m = m
        self.n = len(m)

    def reduce(self):
        (x, y), c = self.unit_edge()
        while c:
            for w in range(self.n):
                for z in range(self.n):
                    if w != x and z != y:
                        self.m[w, z] += -self.m[w, y] * c * self.m[x, z]

            self.m[x, :] = 0
            self.m[y, :] = 0
            self.m[:, x] = 0
            self.m[:, y] = 0

            (x, y), c = self.unit_edge()

    def unit_edge(self):
        for (x, y), c in np.ndenumerate(self.m):
            if c == 1 or c == -1:
                return (x, y), c
        return (False, False), False
