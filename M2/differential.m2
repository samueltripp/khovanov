-- maybe will use for labeling things later
TwoWayTable = new Type of HashTable;
twoWayTable = method();
twoWayTable(HashTable) := (ht) -> (
    out := new TwoWayTable from hashTable {
	(global keyToValue) => new MutableHashTable,
	(global valueToKey) => new MutableHashTable
    };
    
    htKeys := keys(ht);
    for i from 0 to #htKeys - 1 do (
	out.keyToValue#(htKeys#i) = ht#(htKeys#i);
	out.valueToKey#(ht#htKeys#i) = htKeys#i;
    );
    
    return out;
);

-- constructor for chain complexes
DifferentialGradedModule = new Type of HashTable;
differentialGradedModule = method(Options => {squaresToZero => true});
-- g : {Grading}
differentialGradedModule(Module, Matrix, List) := DifferentialGradedModule =>
opts ->
(m, d, g) -> (
    if opts.squaresToZero then assert(d^2 == 0);
    -- TODO: assert that g is compatible with d
    return new DifferentialGradedModule from hashTable {
        (global m) => m,
	(global d) => d,
	(global g) => g
    };
);

-- the identity complex for the tensor product
-- 0 -> R -> 0
withZeroDifferential = method();
withZeroDifferential(Ring) := DifferentialGradedModule =>
(r) -> differentialGradedModule(r^1, map(r^1, r^1, 0), {0});

-- the map v -> (-1)^|v| v
gradedIdentity = method();
gradedIdentity(DifferentialGradedModule) := (dgm) -> (
    return map(dgm.m, dgm.m, (i,j) -> if i == j then (-1)^(lift(dgm.g#i,ZZ)) else 0);
);

-- tensor product of chain complexes
DifferentialGradedModule ** DifferentialGradedModule := DifferentialGradedModule =>
(dgm1, dgm2) -> (
    m := dgm1.m ** dgm2.m;
    d := dgm1.d ** map(dgm2.m, dgm2.m, 1) + gradedIdentity(dgm1) ** dgm2.d;
    g := for i from 0 to #dgm1.g * #dgm2.g - 1 list dgm1.g#(i // #dgm2.g) + dgm2.g#(i % #dgm2.g);
    return differentialGradedModule(m, d, g, squaresToZero => false);
);