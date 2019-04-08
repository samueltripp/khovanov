from C2Minus import *
from pprint import *

# Get a C2Minus object and find the PreFCC
n = 1
crossings = (1,)
m = len(crossings)

C2M = C2Minus(n, crossings)
pfcc = C2M.preFCC(variant=C2Minus.Variant.BAR)

# Define the base ring and relevant ideals and modules

R_def = 'R = QQ[{0}]'.format(','.join(pfcc.rings[0,0,0].cover_ring().variable_names()))
I_defs = {crossing_key: 'I{0} = ideal({1})'.format(str(crossing_key), str(pfcc.rings[crossing_key,0,0].defining_ideal().gens())[1:-1]) for crossing_key in range(2**m)}
M_defs = {crossing_key: 'M{0} = (R^1)/(I{0}*R^1)'.format(str(crossing_key)) for crossing_key in range(2**m)}

# Define the total complex as a module as a sum of modules at each vertex

M_def = ''
v_indices = {}
f_levels = {}
for i, v in enumerate(pfcc.vertices.keys()):
    f_level = pfcc.vertices[v]
    if f_level in f_levels:
        f_levels[f_level].add(v)
    else:
        f_levels[f_level] = {v}
        
    v_indices[v] = i
    crossing_key = v[0]
    if M_def == '':
        M_def = 'M = M{}'.format(str(crossing_key))
    else:
        M_def += ' ++ M{}'.format(str(crossing_key))
		
# Defining the differential for large objects made my emacs crash so I split it up.

# arrow_strings = set()
# for e in pfcc.edges:
#     if isinstance(e.coefficient, sage.rings.quotient_ring_element.QuotientRingElement):
#         e.coefficient = e.coefficient.lift()
#     arrow_strings.add('({1},{0}) => {2}'.format(str(v_indices[e.source]), str(v_indices[e.target]), str(e.coefficient)))    
# del_def = 'del = map(M, M, {{{0}}})'.format(', '.join(arrow_strings))

arrow_strings = set()
for e in pfcc.edges:
    if isinstance(e.coefficient, sage.rings.quotient_ring_element.QuotientRingElement):
        e.coefficient = e.coefficient.lift()
    arrow_strings.add('({1},{0}) => {2}'.format(str(v_indices[e.source]), str(v_indices[e.target]), str(e.coefficient)))    

arrow_strings = list(arrow_strings)
arrows_per_def = 50
del_defs = ['del = map(M, M, {{{0}}})'.format(', '.join(arrow_strings[0:arrows_per_def]))]
next_arrow = arrows_per_def
while next_arrow < len(arrow_strings):
    del_defs.append('del = del + map(M, M, {{{0}}})'.format(', '.join(arrow_strings[next_arrow:next_arrow+arrows_per_def])))
    next_arrow = next_arrow + arrows_per_def
	
# Define the filtration as a list of submodules of M
    
f_component_strings = []

for f_level in range(0, max(f_levels.keys())+2):
    f_component_strings.append('ideal(0)*M_0')
    
for v, f_level in pfcc.vertices.items():
    index = v_indices[v]
    for i in range(0, f_level+1):
        f_component_strings[i] += ' + R*M_{}'.format(index)

C_def = 'C = [{0}]'.format(', '.join(f_component_strings))

# Print everything for copy/paste into Macaulay2

print(R_def)
for k in I_defs:
    print(I_defs[k])
    print(M_defs[k])
print M_def
for d in del_defs:
    print(d)
print C_def