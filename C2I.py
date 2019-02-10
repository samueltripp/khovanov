class C2I: 

    def __init__(self,n,II):
        self.n = n
        self.II = II
        self.polyring = self.initpolyring(n)
        self.ccx = self.initchaincomplex(n,II.permutation(),self.polyring)

    def size(self):
        return self.n
    
    def TypeII(self):
        return self.II
    
    def polyring(self):
        return self.polyring

    def chaincomplex(self):
        return self.ccx


    @staticmethod
    def initpolyring(n):
        nvars = 6*n-4
        var_names = ['U%s'%p for p in range(1,nvars+1)]
        return PolynomialRing(QQ,nvars,var_names)
        

    def initchaincomplex(self,n,bij,R):
        G = AdditiveAbelianGroup([2])
        G0 = G(vector([0]))
        G1 = G(vector([1]))
        gens = R.gens()
        cxlist = []
        revbij = self.revbij(bij)
        
        # compute the first crossing
        m0 = matrix([gens[0]+gens[1]-gens[2*n]-gens[revbij[1]-1]])
        m1 = matrix([gens[0]+gens[1]+gens[2*n]+gens[revbij[1]-1]])
        cxlist.append(ChainComplex({G0:m0,G1:m1}, base_ring = R,grading_group = G, degree = G1, check=False))

        # compute the rest of the top row except the last one
        for i in range(2,n):
            m0 = matrix([gens[2*i-2]+gens[2*i-1]-gens[2*n+2*i-3]-gens[2*n+2*i-2]])
            m1 = matrix([gens[2*i-2]+gens[2*i-1]+gens[2*n+2*i-3]+gens[2*n+2*i-2]])
            C = ChainComplex({G0:m0,G1:m1}, base_ring = R, grading_group = G, degree = G1, check=False)
            cxlist.append(C)

        # compute the end of the top row
        m0 = matrix([gens[2*n-2]+gens[2*n-1]-gens[4*n-3]-gens[6*n-5]])
        m1 = matrix([gens[2*n-2]+gens[2*n-1]+gens[4*n-3]+gens[6*n-5]])
        cxlist.append(ChainComplex({G0:m0,G1:m1},base_ring = R, grading_group = G, degree = G1, check= False))

        # compute the first crossing in the second row
        m0 = matrix([gens[2*n]+gens[2*n+1]-gens[revbij[2]-1]-gens[4*n-2]])
        m1 = matrix([gens[2*n]+gens[2*n+1]+gens[revbij[2]-1]-gens[4*n-2]])
        cxlist.append(ChainComplex({G0:m0,G1:m1},base_ring = R, grading_group=G,degree=G1,check=False))

        # compute the rest of the second row
        for i in range(2,n):
            m0 = matrix([gens[2*n+2*i-2]+gens[2*n+2*i-1]-gens[4*n-2+2*i-2-1]-gens[4*n-2+2*i-2]])
            m1 = matrix([gens[2*n+2*i-2]+gens[2*n+2*i-1]+gens[4*n-2+2*i-2-1]+gens[4*n-2-2*i-2]])
            C = ChainComplex({G0:m0,G1:m1}, base_ring = R, grading_group = G, degree = G1, check=False)
            cxlist.append(C)
            
        # compute the third row
        for i in range(1,n):
            m0 = matrix([gens[4*n-2+2*i-2]+gens[4*n-2+2*i-1]-gens[revbij[2*i+1]-1]-gens[revbij[2*i+2]-1]])
            m1 = matrix([gens[4*n-2+2*i-2]+gens[4*n-2+2*i-1]+gens[revbij[2*i+1]-1]+gens[revbij[2*i+2]-1]])
            C = ChainComplex({G0:m0,G1:m1},base_ring = R, grading_group = G, degree = G1, check = False)
            cxlist.append(C)

        return cxlist

    def revbij(self,bij):
        return {bij[i]:i for i in range(1,2*self.n+1)};
