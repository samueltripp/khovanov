-- Using cohomological conventions, F(0) = C and F(n) = 0 for some n
-- Define C and del first before running.
-- C =
-- del = 
F = p -> C#(min(max(p,0),#C-1))
eta = p -> inducedMap(F(p)/F(p+1), F(p))
A = (r,p) -> kernel inducedMap(F(p)/F(p+r), F(p), del)
Z = (r,p) -> image inducedMap(, A(r,p), eta(p))
B = (r,p) -> image inducedMap(, F(p-r+1), del)
E = (r,p) -> Z(r,p)/B(r,p)
d = (r,p) -> inducedMap(E(r,p+r), E(r,p), del)
-- Can get pages with commands like "trim E(2,0)"