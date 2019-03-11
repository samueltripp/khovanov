from C2Minus import *

C2M = C2Minus(1,(1,))
h = C2M.homology(5)
print(h.vertices)

C2M = C2Minus(2,(2,2,2))
h = C2M.homology(3)
print(h.vertices)
