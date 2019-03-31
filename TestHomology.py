from C2Minus import *

C2M = C2Minus(2,(2,2,2))
pfcc = C2M.preFCC(variant = C2Minus.Variant.BAR)
ideal_list = list({ring_defining_ideal(ring) for ring in pfcc.rings.values()})
quotient_rings = [pfcc.R.quotient(ideal) for ideal in ideal_list]
a = gen_list_helper(4,(3,quotient_rings[3]))
