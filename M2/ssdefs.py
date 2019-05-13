from C2Minus import *
from pprint import *

n = 2
crossings = (-2,-2,-2)
m = len(crossings)

print('Computing C2Minus...')
C2M = C2Minus(n, crossings)
pfcc = C2M.preFCC(variant=C2Minus.Variant.BAR)
print('PreFCC')

R_def = 'R = QQ[{0}]'.format(','.join(pfcc.rings[0,0,0].cover_ring().variable_names()))
I_defs = {crossing_key: 'I{0} = ideal({1})'.format(str(crossing_key), str(pfcc.rings[crossing_key,0,0].defining_ideal().gens())[1:-1]) for crossing_key in range(2**m)}
M_defs = {crossing_key: 'M{0} = (R^1)/(I{0}*R^1)'.format(str(crossing_key)) for crossing_key in range(2**m)}

M_def = ''
v_indices = {}
f_levels = {}
for i, v in enumerate(sorted(pfcc.vertices.keys())):
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
    
Gr_defs = ['Gr = new MutableList']

for f_level in range(0, max(f_levels.keys())+1):
    Gr_defs.append('Gr#{} = 0*M'.format(f_level))
    
for v, f_level in pfcc.vertices.items():
    Gr_defs[f_level+1] += ' + R*M_{}'.format(v_indices[v])

print(R_def)
for k in I_defs:
    print(I_defs[k])
    print(M_defs[k])
print M_def
for k in Gr_defs:
    print k
for d in del_defs:
    print(d)