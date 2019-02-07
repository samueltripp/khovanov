class Braid:
    # n - the number of strands in the braid
    # word - a word in the braid group B_n, given as a tuple {-(n-1),...,-1,1,...,(n-1)}*
    #
    # conventions: braids are pictured vertically, strands are numbered 1-n left-to-right,
    # and the generator i (-i) moves the ith strand under (over) the (i+1)st strand.
    def __init__(self, n, word):
        for i in word:
            assert 1 <= abs(i) < n, "Invalid word"

        self.n = n
        self.word = word

    # returns a dict representing image of this braid under the canonical projection B_n -> S_n
    def permutation(self):
        p = {i: i for i in range(1, self.n+1)}
        if (self.word == {}): return p
        for swap in reversed(self.word):
            swap = abs(swap)
            p[swap], p[swap+1] = p[swap+1], p[swap]
        return p
